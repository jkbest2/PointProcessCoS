library(INLA)          # Used for model fitting, mesh generation etc.
library(inlabru)       # More convenient interface to INLA, allows for nonlinear
                       # effects
library(spatstat)      # Provides rLGCP
library(RandomFields)  # Defines Matern correlation functions
library(rgeos)         # Useful for polygon intersections etc.
## This file provides functions for creating the dual mesh as well as other
## conveniences. It is taken from the .zip file available here:
## (http://www.r-inla.org/spde-book)
source("spde-book/spde-book-functions.R")
source("ppcos_functions.R")
## Load tidyverse last so that purrr::map isn't masked by maps::map
library(tidyverse)     # Data management etc.
library(gridExtra)
library(viridis)
library(TMB)

set.seed(1221)

## Intermediate plots?
inter_plots <- FALSE

###=============================================================================
### Simulate data --------------------------------------------------------------
## Use a window that is 100×100. Create an owin object for the `rLGCP` function,
## and an `sp::SpatialPolygons` object for masking later.
xrange <- yrange <- c(0, 100)
ppwin <- owin(xrange = xrange, yrange = yrange)
pp_domain <- matrix(c(xrange[1], yrange[1],
                      xrange[2], yrange[1],
                      xrange[2], yrange[2],
                      xrange[1], yrange[2],
                      xrange[1], yrange[1]),
                    ncol = 2L, byrow = TRUE)
pp_dom_poly <- SpatialPolygons(list(Polygons(list(Polygon(pp_domain)), "0")))

## Generate a log-Gaussian Cox process realization. Save the intensity surface.
## This intensity surface is what we are ultimately trying to make inference on.
rho = 50
sigma2 = 1 ^ 2
## Convert to the Matérn parameterization used by INLA
kappa2 <- 8 / rho ^ 2
tau <- sqrt(4 * pi * kappa2 * sigma2)
## Want to make the underlying mean high enough that we have lots of points.
## Otherwise our point process data won't be very informative. For mu = -2, we
## expect to see exp(-3) * 100^2 = 498 points.
mu = -2
## Set nu = 1; this corresponds to alpha = 2 in the SPDE formulation.
pp1 <- rLGCP(model = "matern", nu = 1,
             mu = mu, var = sigma2, scale = rho,
             win = ppwin, eps = 1, saveLambda = TRUE)
## Extract the point process locations as different objects
pp_sp <- data.frame(pp1)
coordinates(pp_sp) <- c("x", "y")
## Total number of points
n_pp <- pp1$n

### Simulate the area-wide estimate
## area_cv <- 0.25
## area_sd <- sqrt(log(area_cv ^ 2 + 1))
## area_est <- rlnorm(1, log(sum(attr(pp1, "Lambda")$v)), area_sd)
area_est <- rpois(1, sum(attr(pp1, "Lambda")$v))

###=============================================================================
### Prepare data and mesh ------------------------------------------------------

### This follows the outline of an example set forth in *Advanced Spatial
### Modeling with Stochastic Partial Differential Equations Using R and INLA*
### Chapter 4, available here:
### (https://becarioprecario.bitbucket.io/spde-gitbook/ch-lcox.html). It follows
### Simpson 2016 (https://doi.org/10.1093/biomet/asv064) to construct an
### efficient representation of the log-Gaussian Cox process by using the dual
### of the mesh that is constructed when using the SPDE representation of the
### Matern spatial field.

## Generate mesh for FEM representation of SPDE. The `loc` argument gives a
## starting point for the vertex locations, `loc.domain` defines the domain of
## the point process, `max.edge` controls the maximum edge lengths in the mesh
## (so smaller values result in a finer mesh at the expense of more vertices to
## fit), `cutoff` controls how close two `loc`s are before they are merged to a
## single vertex, and `offset` controls how far the inner mesh extends from
## `loc.domain` and the outer mesh from there. The `loc` argument is optional if
## the others are specified.
mesh <- inla.mesh.2d(loc.domain = pp_domain,
                     max.edge = c(5, 20),
                     cutoff = 1,
                     offset = c(5, 20))
## We'll need the number of vertices in the mesh later for fitting.
n_vert <- mesh$n

if (inter_plots) {
  ## Visually check the mesh. Triangles should not be "long and skinny", the finer
  ## inner mesh should extend beyond the green square, and the outer boundary of
  ## the mesh should extend substantially further from there to compensate for
  ## edge effects (inflated variance at the edges).
  ggplot() +
    gg(mesh) +
    gg(pp_dom_poly) +
    gg(pp_sp, alpha = 0.5) +
    coord_fixed()
}

## Because we are estimating using maximum likelihood here, it doesn't really
## matter whether we use `inla.spde2.pcmatern` or `inla.spde2.matern` to
## construct the SPDE object that we use to get a precision matrix.
spde <- inla.spde2.pcmatern(mesh = mesh,
                            prior.range = c(20, 0.01),
                            prior.sigma = c(5, 0.01))
## Then specify the precision matrix using the generative parameters from above.
## This returns a sparse matrix that we can pass to TMB.
spde_Q <- inla.spde2.precision(spde,
                               theta = c(log(rho), log(sqrt(sigma2))))

## Create the dual mesh (divide the domain into polygons according to the
## closest vertex). This uses the `book.mesh.dual` function `source`d in above.
mesh_dual <- book.mesh.dual(mesh)

### Get quadrat partially-observed point process observations sorted
## Number of quadrats to map the point process
## n_blocks <- 5L
## List of the lower-left corner coordinates for the blocks, to generate the
## coordinates of all possible blocks.
block_ll <- seq(0, 75, 25)
blocks <- cross2(block_ll, block_ll) %>%
  map(lift(c)) %>%
  map(block_poly)

make_poly <- function(x) {
  Polygon(matrix(x, ncol = 2, byrow = TRUE))
}

polys <- list(make_poly(c(10, 85,
                          25, 90,
                          35, 80,
                          50, 70,
                          30, 65,
                          20, 70,
                          10, 85)),
              make_poly(c(75, 15,
                          60, 15,
                          55, 30,
                          60, 40,
                          75, 50,
                          80, 40,
                          70, 30,
                          75, 15)),
              make_poly(c(65, 85,
                          85, 85,
                          85, 65,
                          65, 65,
                          65, 85)),
              make_poly(c(15, 15,
                          10, 20,
                          30, 40,
                          35, 35,
                          15, 15)))

## Choose a sample of quadrats from within the coarse 16×16 grid, then get the
## dual mesh areas within each quadrat and the point counts within each dual
## mesh/quadrat intersection. These are used to calculate the likelihood of the
## partially observed point process.
## quad_poly <- blocks[c(10)]
## quad_poly <- append(quad_poly, irr_poly)

## sp <- SpatialPolygons(c(list(Polygons(irr_poly, "irr_quads")),
##                         list(Polygons(quad_poly, "reg_quads"))))
poly_sp <- SpatialPolygons(list(Polygons(polys, "polys")))
ggplot() +
  gg(mesh) +
  gg(pp_dom_poly) +
  gg(poly_sp)

## quadrat_sp <- SpatialPolygons(list(Polygons(quad_poly, "quadrats")))
quadrat_sp <- poly_sp
quad_wts <- get_wts(quadrat_sp, mesh_dual)
quad_counts <- count_points(pp_sp, mesh_dual, quadrat_sp)

## We also need the area of each dual mesh polygon in order to do the numerical
## integration to get the total intensity for comparison with the area-wide
## abundance estimate.
dual_wts <- get_wts(pp_dom_poly, mesh_dual)

###=============================================================================
### Fit the model --------------------------------------------------------------

## pmarg_range <- qpois(c(0.05, 0.995), sum(attr(pp1, "Lambda")$v))
## pmarg_vec <- pmarg_range[1]:pmarg_range[2]
## N_pmarg <- length(pmarg_vec)

data <- list(N_vertices = mesh$n,
             quadrat_count = quad_counts,
             quadrat_exposure = quad_wts,
             dual_exposure = dual_wts,
             area_est = area_est,
             ## area_sd = area_sd,
             ## N_pmarg = N_pmarg,
             ## pmarg_range = pmarg_vec,
             Q = spde_Q)
pars <- list(mu = mu,
             spat = rep(0, mesh$n))

compile("ppcos_quads.cpp")
dyn.load(dynlib("ppcos_quads"))

obj <- MakeADFun(data, pars,
                 random = "spat",
                 DLL = "ppcos_quads")

fit <- optim(obj$par, obj$fn, obj$gr, method = "BFGS",
             control = list(maxit = 1000L))

h <- optimHess(fit$par, obj$fn, obj$gr)

rep <- obj$report()
sdr <- sdreport(obj)

###=============================================================================
### Figures etc. ---------------------------------------------------------------

pp_sp$obs <- ifelse(is.na(over(pp_sp, quadrat_sp)),
                    FALSE, TRUE)

pxl <- pixels(mesh, nx = 150, ny = 150, mask = pp_dom_poly)
pxl_A <- inla.spde.make.A(mesh, pxl)
lin_est <- (rep$mu + pxl_A %*% rep$spat)[, 1]
## To get the standard deviation of the linear predictor at each pixel location,
## we take the diagonal of A Sigma A', add the variance of the mean, and then
## take the diagonal. This gives the variance of each pixel, which we take the
## square root of.
lin_sd <- sqrt(diag(pxl_A %*% Diagonal(x = sdr$diag.cov.random) %*% t(pxl_A)) +
               sdr$cov.fixed[1, 1])
lam_est <- exp(lin_est)

res_df <- tibble(x = coordinates(pxl)[, 1],
                 y = coordinates(pxl)[, 2],
                 lin_est = lin_est,
                 lin_sd = lin_sd,
                 lam_est = lam_est,
                 lam_true = c(t(attr(pp1, "Lambda")$v)),
                 lam_diff = lam_est - lam_true,
                 lin_diff = lin_est - log(lam_true))

## Left plot, true surface with observations
true_plot <- ggplot(res_df, aes(x = x, y = y)) +
  geom_tile(aes(fill = lam_true)) +
  scale_fill_viridis_c(option = "D", limits = c(0.01, 0.91)) +
  coord_fixed() +
  geom_point(data = as.data.frame(pp_sp[pp_sp$obs, ]),
             alpha = 0.5, size = 0.2) +
  geom_path(data = fortify(quadrat_sp),
            aes(x = long, y = lat, group = group),
            alpha = 0.5) +
  labs(fill = expression(Lambda(s)),
       title = "(a)") +
  guides(fill = FALSE) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank())

## Right plot, reconstructed surface
est_plot <- ggplot(res_df, aes(x = x, y = y)) +
  geom_tile(aes(fill = lam_est)) +
  scale_fill_viridis_c(option = "D", limits = c(0.01, 0.91)) +
  coord_fixed() +
  ## geom_point(data = as.data.frame(pp_sp[pp_sp$obs, ]),
  ##            alpha = 0.5) +
  geom_path(data = fortify(quadrat_sp),
            aes(x = long, y = lat, group = group, fill = NULL),
            alpha = 0.5) +
  labs(fill = expression(Lambda(s)),
       title = "(b)") +
  guides(fill = FALSE) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank())

combo_plot <- grid.arrange(true_plot, est_plot, ncol = 2)

ggsave("ppcos_plot.pdf", combo_plot, width = 6, height = 4)

## Alternate plots -------------------------------------------------------------
ggplot(res_df, aes(x = x, y = y)) +
  geom_tile(aes(fill = lam_est)) +
  scale_fill_viridis_c(option = "D") +
  coord_fixed() +
  geom_contour(aes(z = lam_true),
               breaks = seq(0, 1, 0.25)) +
  geom_path(data = fortify(quadrat_sp),
            aes(x = long, y = lat, group = group, fill = NULL))

ggplot(res_df, aes(x = x, y = y, fill = lam_est)) +
  geom_tile(alpha = 0.5) +
  scale_fill_viridis_c(option = "D") +
  coord_fixed() +
  geom_point(data = as.data.frame(pp_sp),
             aes(alpha = obs, fill = NULL), guide = FALSE) +
  geom_path(data = fortify(quadrat_sp),
            aes(x = long, y = lat, group = group, fill = NULL))

## Plot the standard deviation of the linear predictor (log intensity) estimate.
## It's noisy, probably because of the discretization, but you can see a
## reduction in standard deviation around the observed quadrats.
ggplot(res_df, aes(x = x, y = y, fill = lin_sd)) +
  geom_tile() +
  scale_fill_viridis_c(option = "E") +
  coord_fixed() +
  geom_path(data = fortify(quadrat_sp),
            aes(x = long, y = lat, group = group, fill = NULL))

ggplot(res_df, aes(x = x, y = y)) +
  geom_tile(aes(fill = lin_diff)) +
  ## scale_fill_viridis_c( option = "E") +
  scale_fill_distiller(type = "div",
                       limits = c(-2.5, 2.5)) +
  coord_fixed() +
  geom_path(data = fortify(quadrat_sp),
            aes(x = long, y = lat, group = group)) +
  geom_contour(aes(z = lin_diff), breaks = seq(-2, 2, 0.5))

