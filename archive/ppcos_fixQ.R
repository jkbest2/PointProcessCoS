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
library(TMB)

set.seed(1221)

### ============================================================================
### Define domain of interest
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

###=============================================================================
### Simulate point process -----------------------------------------------------
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
mu = -3
## Set nu = 1; this corresponds to alpha = 2 in the SPDE formulation.
pp1 <- rLGCP(model = "matern", nu = 1,
             mu = mu, var = sigma2, scale = rho,
             win = ppwin, eps = 1, saveLambda = TRUE)
## Extract the point process locations as different objects
pp_df <- data.frame(pp1)
coordinates(pp_df) <- c("x", "y")
## Total number of points
n_pp <- pp1$n

### ============================================================================
### Prepare data and mesh-------------------------------------------------------
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

## Visually check the mesh. Triangles should not be "long and skinny", the finer
## inner mesh should extend beyond the green square, and the outer boundary of
## the mesh should extend substantially further from there to compensate for
## edge effects (inflated variance at the edges).
ggplot() +
  gg(mesh) +
  gg(pp_dom_poly) +
  gg(pp_df, alpha = 0.5) +
  coord_fixed()

## Because we are estimating using maximum likelihood here, it doesn't really
## matter whether we use `inla.spde2.pcmatern` or `inla.spde2.matern` to
## construct the SPDE object that we use to get a precision matrix.
spde <- inla.spde2.pcmatern(mesh = mesh,
                            prior.range = c(20, 0.01),
                            prior.sigma = c(5, 0.01))
## Then specify the precision matrix using the generative parameters from above.
## This returns a sparse matrix that we can pass to TMB directly.
spde_Q <- inla.spde2.precision(spde,
                               theta = c(log(rho), log(sqrt(sigma2))))

## Create the dual mesh (divide the domain into polygons according to the
## closest vertex). This uses the `book.mesh.dual` function `source`d in above.
mesh_dual <- book.mesh.dual(mesh)

### ============================================================================
### Simulate areal observations ------------------------------------------------
## Observations of total area population. Currently just independent Poisson
## draws from the total intensity; not clear exactly where this would come from
## in a realistic scenario.
n_areaobs <- 3L
## Each year's total number of nests is drawn from a Poisson distribution with
## rate equal to the total intensity of the domain
tot_intens <- sum(attr(pp1, "Lambda")$v)
Nt <- rpois(n_areaobs, tot_intens)
## Observations are log-normal, with a CV of 0.25
Nhat <- rlnorm(n_areaobs, log(Nt), sqrt(log(0.25^2 + 1)))

## Quick plot of the point process and underlying intensity surface
## Matrix::image(attr(pp1, "Lambda"))
## points(pp1)

### ============================================================================
### Delineate blocks for point-process observations ----------------------------
## Number of blocks to simulate observations
n_blocks <- 5L
## List of the lower-left corner coordinates for the blocks, to generate the
## coordinates of all possible blocks.
block_ll <- seq(0, 75, 25)
blocks <- cross2(block_ll, block_ll) %>%
  map(lift(c)) %>%
  map(block_poly)

## Get the total intensity for each observed block, then simulate a
## Poisson-distributed observation for each block based on the total intensity
## in each block.
block_df <- tibble(
  Poly = blocks[sample(1:16, n_blocks, replace = FALSE)]) %>%
  mutate(intensity = map_dbl(Poly, ~ total_intens(., attr(pp1, "Lambda"))),
         obs_count = map_dbl(intensity, ~ rpois(1, .))),
         pp_wts = map(Poly, get_wts, mesh_dual = mesh_dual),
         pp_counts = map(Poly, function(poly) count_points(pp1, mask = poly)))

quadrat_sp <- SpatialPolygons(list(Polygons(block_df$Poly, "quadrats")))
quad_wts <- get_wts(quadrat_sp, mesh_dual)


## Get the areas withing the domain for each polygon in the dual mesh
wts <- get_wts(pp_dom_poly, mesh_dual)
## Check to make sure that the total area calculated above is equal to the known
## total area of the point process domain.
if (sum(wts) != diff(xrange) * diff(yrange)) stop("Weights wrong!")

## Get the point counts for each dual mesh cell
pp_sp <- SpatialPoints(coords = cbind(pp1$x, pp1$y))

cell_counts <- sapply(1:length(mesh_dual), function(i) {
  length(gIntersection(mesh_dual[i, ], pp_sp))
})
## Check to make sure sum of cell counts matches total count in point process
## realization
if (sum(cell_counts) != n_pp) stop("Point counts wrong!")

## Plot the dual mesh, marking which polygons are included in the point process
## domain. The `wts` vector above contains the area of each polygon that is
## within the domain (blue box).
## plot(mesh$loc, asp = 1, col = (wts == 0) + 1, pch = 19, xlab = "", ylab = "")
## plot(mesh_dual, add = TRUE)
## lines(pp_domain, col = "blue", lwd = 2)

data <- list(N_vertices = mesh$n,
             N_areaobs = n_areaobs,
             pp_counts = cell_counts,
             pp_cell_area = wts,
             areaobs = areaobs,
             Q = spde_Q)
pars <- list(mu = mu,
             spat = rep(0, mesh$n))

compile("ppcos_fixQ.cpp")
dyn.load(dynlib("ppcos_fixQ"))

obj <- MakeADFun(data, pars,
                 DLL = "ppcos_fixQ")

fit <- optim(obj$par, obj$fn, obj$gr, method = "BFGS",
             control = list(maxit = 1000L))

h <- optimHess(fit$par, obj$fn, obj$gr)

rep <- obj$report()
sdr <- sdreport(obj)

