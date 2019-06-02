library(INLA)          # Used for model fitting, mesh generation etc.
library(inlabru)       # More convenient interface to INLA, allows for nonlinear
                       # effects
library(spatstat)      # Provides rLGCP
library(RandomFields)  # Defines Matern correlation functions
library(rgeos)         # Useful for polygon intersections etc.
library(maptools)
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
inter_plots <- TRUE

### ============================================================================
### Define domain --------------------------------------------------------------
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

### Define mesh ----------------------------------------------------------------
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
    coord_fixed()
}

## Define pixels for plotting and a projection matrix to those pixels
pxl <- pixels(mesh, nx = 150, ny = 150, mask = pp_dom_poly)
pxl_A <- inla.spde.make.A(mesh, pxl)
pxl_df <- tibble(x = coordinates(pxl)[, 1],
                 y = coordinates(pxl)[, 2])

## Create the dual mesh (divide the domain into polygons according to the
## closest vertex). This uses the `book.mesh.dual` function `source`d in above.
mesh_dual <- book.mesh.dual(mesh)
if (inter_plots) {
  ggplot() +
    gg(mesh_dual) +
    gg(pp_dom_poly)
}

### ============================================================================
### Generate a spatially varying covariate -------------------------------------
cov_rho = 75
cov_sigma2 = 1 ^ 2
spde <- inla.spde2.pcmatern(mesh = mesh,
                            prior.range = c(20, 0.01),
                            prior.sigma = c(5, 0.01))
cov_Q <- inla.spde2.precision(spde,
                              theta = c(log(cov_rho),
                                        log(sqrt(cov_sigma2))))
spat_cov <- backsolve(chol(cov_Q), rnorm(n_vert))

## Create design matrix for covariate effect (at each vertex location)
X <- cbind(rep(1, n_vert), spat_cov)

if (inter_plots) {
  pxl_df %>%
    mutate(covar = (pxl_A %*% spat_cov)[, 1]) %>%
    ggplot(aes(x = x, y = y, fill = covar)) +
    geom_tile() +
    scale_fill_viridis_c(option = "E") +
    coord_fixed()
}

## Calculate true mean vector
mu_pxl <- (-3.0 + 1 * pxl_A %*% spat_cov)[, 1]
## Convert to "im" object for `rLGCP` below
mu_im <- as.im(matrix(mu_pxl,
                      nrow = 100, ncol = 100,
                      byrow = TRUE),
               W = ppwin)
if (inter_plots) {
  pxl_df %>%
    mutate(mu = mu_pxl) %>%
    ggplot(aes(x = x, y = y, fill = mu)) +
    geom_tile() +
    scale_fill_viridis_c(option = "E") +
    coord_fixed()
}

### ============================================================================
### Generate point process realization -----------------------------------------
## Generate a log-Gaussian Cox process realization. Save the intensity surface.
## This intensity surface is what we are ultimately trying to make inference on.
## Here I set the range to be much larger than the domain, and the variance to
## be very small to essentially remove the spatial component, leaving just the
## mean structure. This simplifies the LGCP to a Cox process to match the rest
## of the text more closely
rho <- 500
sigma2 <- 1e-10 ^2
spde_Q <- inla.spde2.precision(spde, theta = c(log(rho), log(sqrt(sigma2))))
## Convert to the Matérn parameterization used by INLA
kappa2 <- 8 / rho ^ 2
tau <- sqrt(4 * pi * kappa2 * sigma2)
## Want to make the underlying mean high enough that we have lots of points.
## Otherwise our point process data won't be very informative. For mu = -2, we
## expect to see exp(-3) * 100^2 = 498 points.
## mu <- -2
## mu <- function(x, y) -3.25 + 0.02 * y + 0.02 * x
## Set nu = 1; this corresponds to alpha = 2 in the SPDE formulation.
pp1 <- rLGCP(model = "matern", nu = 1,
             mu = mu_im, var = sigma2, scale = rho,
             win = ppwin, eps = 1, saveLambda = TRUE)
## Extract the point process locations as different objects
pp_sp <- data.frame(pp1)
coordinates(pp_sp) <- c("x", "y")
## Total number of points
n_pp <- pp1$n

if (inter_plots) {
  image(attr(pp1, "Lambda"))
  points(pp1)
}

### Simulate the areal estimates. Define polygons to use, get the total
### intensity within that polygon, then simulates Poisson observations.
area_poly_list <- list(
  make_poly(c(0, 0,
              0, 50,
              50, 50,
              50, 0,
              0, 0)),
  make_poly(c(50, 0,
              50, 50,
              100, 50,
              100, 0,
              50, 0)),
  make_poly(c(0, 50,
              0, 100,
              50, 100,
              50, 50,
              0, 50)),
  make_poly(c(50, 50,
              50, 100,
              100, 100,
              100, 50,
              50, 50)))
area_polys_list <- lapply(1:length(area_poly_list),
                          function(idx) Polygons(list(area_poly_list[[idx]]),
                                                paste0("area", idx)))
area_poly <- SpatialPolygons(area_polys_list)

intens_sgdf <- as.SpatialGridDataFrame.im(attr(pp1, "Lambda"))
area_intens <- over(area_poly, intens_sgdf, fn = sum)
area_est <- rpois(4, area_intens$v)
area_wts <- vapply(area_polys_list,
                   function(poly) {
                     get_wts(SpatialPolygons(list(poly)), mesh_dual = mesh_dual)
                   }, FUN.VALUE = rep(0.0, n_vert))

### Get partially-observed point process observations sorted. Use arbitrary
### polygons to illustrate flexibility of approach.
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

poly_sp <- SpatialPolygons(list(Polygons(polys, "polys")))
if (inter_plots) {
  ggplot() +
    gg(mesh) +
    gg(pp_dom_poly) +
    gg(poly_sp)
}

quadrat_sp <- poly_sp
quad_wts <- get_wts(quadrat_sp, mesh_dual)
quad_counts <- count_points(pp_sp, mesh_dual, quadrat_sp)

## We also need the area of each dual mesh polygon in order to do the numerical
## integration to get the total intensity for comparison with the area-wide
## abundance estimate.
dual_wts <- get_wts(pp_dom_poly, mesh_dual)

###=============================================================================
### Fit the model --------------------------------------------------------------
data <- list(N_vertices = mesh$n,
             N_areal = length(area_poly_list),
             X = X,
             quadrat_count = quad_counts,
             quadrat_exposure = quad_wts,
             area_exposure = area_wts,
             area_est = area_est,
             exposure = dual_wts,
             area_est = area_est)
pars <- list(beta = c(-3.0, 1.0))

compile("ppcos_quads.cpp")
dyn.load(dynlib("ppcos_quads"))

obj <- MakeADFun(data, pars,
                 DLL = "ppcos_quads")

fit <- optim(obj$par, obj$fn, obj$gr, method = "BFGS",
             control = list(maxit = 1000L))

h <- optimHess(fit$par, obj$fn, obj$gr)

rep <- obj$report()
sdr <- sdreport(obj)

### Fit the model with no pp obs ------------------------------------------------
dat2 <- list(N_vertices = mesh$n,
             N_areal = length(area_poly_list),
             X = X,
             quadrat_count = quad_counts,
             quadrat_exposure = rep(0.0, mesh$n),
             area_exposure = area_wts,
             area_est = area_est,
             exposure = dual_wts,
             area_est = area_est)
pars <- list(beta = c(-3.0, 1.0))

obj2 <- MakeADFun(dat2, pars,
                  DLL = "ppcos_quads")

fit2 <- optim(obj2$par, obj2$fn, obj2$gr, method = "BFGS",
              control = list(maxit = 1000L))

h2 <- optimHess(fit2$par, obj2$fn, obj2$gr)

rep2 <- obj2$report()
sdr2 <- sdreport(obj2)

###=============================================================================
### Figures etc. ---------------------------------------------------------------

pp_sp$obs <- ifelse(is.na(over(pp_sp, quadrat_sp)),
                    FALSE, TRUE)

pxl <- pixels(mesh, nx = 150, ny = 150, mask = pp_dom_poly)
pxl_A <- inla.spde.make.A(mesh, pxl)
lin_est <- (pxl_A %*% (rep$mu))[, 1]
## To get the standard deviation of the linear predictor at each pixel location,
## we take the diagonal of A Sigma A', add the variance of the mean, and then
## take the diagonal. This gives the variance of each pixel, which we take the
## square root of.
## lin_sd <- sqrt(diag(pxl_A %*% Diagonal(x = sdr$diag.cov.random) %*% t(pxl_A)) +
##                sdr$cov.fixed[1, 1])
lam_est <- exp(lin_est)

res_df <- tibble(x = coordinates(pxl)[, 1],
                 y = coordinates(pxl)[, 2],
                 covar = mu_pxl,
                 lin_est = lin_est,
                 ## lin_sd = lin_sd,
                 lam_est = lam_est,
                 lam_true = c(t(attr(pp1, "Lambda")$v)),
                 lam_diff = lam_est - lam_true,
                 lin_diff = lin_est - log(lam_true))

covar_plot <- ggplot(res_df, aes(x = x, y = y)) +
  geom_tile(aes(fill = covar)) +
  scale_fill_viridis(option = "A") +
  labs(title = "(a)") +
  guides(fill = FALSE) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank())

## Get min and max values of lambda, both true and fitted to use to set limits
## on the color scale
col_lims <- c(min(c(res_df$lam_est, res_df$lam_true)),
              max(c(res_df$lam_est, res_df$lam_true)))

## Left plot, true surface with observations
true_plot <- ggplot(res_df, aes(x = x, y = y)) +
  geom_tile(aes(fill = lam_true)) +
  scale_fill_viridis_c(option = "D", limits = col_lims) +
  coord_fixed() +
  geom_point(data = as.data.frame(pp_sp[pp_sp$obs, ]),
             alpha = 0.5, size = 0.2) +
  geom_path(data = fortify(quadrat_sp),
            aes(x = long, y = lat, group = group),
            alpha = 0.5) +
  labs(fill = expression(Lambda(s)),
       title = "(b)") +
  guides(fill = FALSE) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank())

## Right plot, reconstructed surface
est_plot <- ggplot(res_df, aes(x = x, y = y)) +
  geom_tile(aes(fill = lam_est)) +
  scale_fill_viridis_c(option = "D", limits = col_lims) +
  coord_fixed() +
  ## geom_point(data = as.data.frame(pp_sp[pp_sp$obs, ]),
  ##            alpha = 0.5) +
  geom_path(data = fortify(quadrat_sp),
            aes(x = long, y = lat, group = group, fill = NULL),
            alpha = 0.5) +
  labs(fill = expression(Lambda(s)),
       title = "(c)") +
  guides(fill = FALSE) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank())

combo_plot <- grid.arrange(covar_plot, true_plot, est_plot,
                           nrow = 2, ncol = 2)

ggsave("ppcos_plot.png", combo_plot, width = 6, height = 6)

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

## ggplot(res_df, aes(x = x, y = y)) +
##   geom_tile(aes(fill = lin_diff)) +
##   ## scale_fill_viridis_c( option = "E") +
##   scale_fill_distiller(type = "div",
##                        limits = c(-2.5, 2.5)) +
##   coord_fixed() +
##   geom_path(data = fortify(quadrat_sp),
##             aes(x = long, y = lat, group = group)) +
##   geom_contour(aes(z = lin_diff), breaks = seq(-2, 2, 0.5))

