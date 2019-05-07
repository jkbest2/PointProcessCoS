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
## Load tidyverse last so that purrr::map isn't masked by maps::map
library(tidyverse)     # Data management etc.

set.seed(1111)

###=============================================================================
### Simulate data --------------------------------------------------------------
## These ranges are set as `c(0, 128)` because that is the discretization that
## `rLGCP` uses for the Lambda matrix. This way each cell in the lambda matrix
## is area one and comparison with our results later is slightly more
## convenient.
xrange <- yrange <- c(0, 128)
## yrange <- c(0, 128)
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
## Want to make the underlying mean high enough that we have lots of points.
## Otherwise our point process data won't be very informative. For mu = -2, we
## expect to see exp(-2) * 128^2 = 2217 points.
mu = -2
## Set nu = 1; this corresponds to alpha = 2 in the SPDE formulation.
pp1 <- rLGCP(model = "matern", nu = 1,
             mu = mu, var = sigma2, scale = rho,
             win = ppwin, eps = 1, saveLambda = TRUE)
## Extract the point process locations as different objects
pp_df <- data.frame(pp1)
coordinates(pp_df) <- c("x", "y")
pp_locs <- data.frame(s1 = pp1$x, s2 = pp1$y)
pp_mat <- as.matrix(pp_locs)
n_pp <- pp1$n

## Quick plot of the point process and underlying intensity surface
Matrix::image(attr(pp1, "Lambda"))
points(pp1)

## Function to generate a Polygon of a block from the lower-left corner
block_poly <- function(ll, w = 32) {
  coords <- matrix(c(ll[1],     ll[2],
                     ll[1] + w, ll[2],
                     ll[1] + w, ll[2] + w,
                     ll[1],     ll[2] + w,
                     ll[1],     ll[2]),
                   byrow = TRUE, ncol = 2L)
  Polygon(coords, hole = FALSE)
}
## Number of blocks to simulate observations
n_blocks <- 5L
## List of the lower-left corner coordinates for the blocks, to generate the
## coordinates of all possible blocks.
block_ll <- seq(0, 96, 32)
blocks <- cross2(block_ll, block_ll) %>%
  map(lift(c)) %>%
  map(block_poly)

## Get elements of intensity that are within a block. Takes a Polygon and the
## saved "Lambda" attribute of the log-Gaussian Cox process generated above.
## This is a rough integral over each block. The `eps` argument in the `rLGCP`
## call above controls the grid spacing and so the integral approximation.
total_intens <- function(blk_poly, Lambda) {
  crds <- blk_poly@coords
  x_idx <- which((Lambda$xcol > crds[1, 1]) & (Lambda$xcol < crds[3, 1]))
  y_idx <- which((Lambda$yrow > crds[1, 2]) & (Lambda$yrow < crds[3, 2]))
  sum(Lambda$v[x_idx, y_idx])
}

## Get the total intensity for each observed block, then simulate a
## Poisson-distributed observation for each block based on the total intensity
## in each block.
block_df <- tibble(
  Poly = blocks[sample(1:16, n_blocks, replace = FALSE)]) %>%
  mutate(intensity = map_dbl(Poly, ~ total_intens(., attr(pp1, "Lambda"))),
         obs_count = map_dbl(intensity, ~ rpois(1, .)))

###=============================================================================
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
mesh <- inla.mesh.2d(#loc = pp_loc_mat,
                     loc.domain = pp_domain,
                     max.edge = c(8, 20),
                     cutoff = 1,
                     offset = c(5, 20))
## We'll need the number of vertices in the mesh later for fitting.
n_vert <- mesh$n

## Plot the resulting mesh. Triangles should not be "long and skinny", the finer
## inner mesh should extend beyond the green square, and the outer boundary of
## the mesh should extend substantially further from there to compensate for
## edge effects (inflated variance at the edges).
ggplot() +
  gg(mesh) +
  gg(pp_dom_poly) +
  gg(pp_df, alpha = 0.5) +
  coord_fixed()

## Declare a spatial field using the mesh, placing penalized complexity priors
## (Fuglstad 2018, https://doi.org/10.1080/01621459.2017.1415907) on the spatial
## range and standard deviation. These declare that Pr(range < 10) = 0.01 and
## Pr(sigma > 5) = 0.01, placing substantial mass on the area around our
## generative values without being overly restrictive.
spde <- inla.spde2.pcmatern(mesh = mesh,
                            prior.range = c(20, 0.01),
                            prior.sigma = c(5, 0.01))

## Create the dual mesh (divide the domain into polygons according to the
## closest vertex). This uses the `book.mesh.dual` function `source`d in above.
mesh_dual <- book.mesh.dual(mesh)

## Calculate the weights (areas) of each polygon of the dual mesh using `rgeos`
## functions, intersected with the point process domain. Taken from Chapter 4 of
## Advanced Spatial Modeling (linked above).
wts <- sapply(1:length(mesh_dual), function(i) {
  if (gIntersects(mesh_dual[i, ], pp_dom_poly)) {
    area <- gArea(gIntersection(mesh_dual[i, ], pp_dom_poly))
  } else {
    area <- 0
  }
  area
})
## Check to make sure that the total area calculated above is equal to the known
## total area of the point process domain.
if (sum(wts) != diff(xrange) * diff(yrange)) stop("Weights wrong!")

## Plot the dual mesh, marking which polygons are included in the point process
## domain. The `wts` vector above contains the area of each polygon that is
## within the domain (blue box).
plot(mesh$loc, asp = 1, col = (wts == 0) + 1, pch = 19, xlab = "", ylab = "")
plot(mesh_dual, add = TRUE)
lines(pp_domain, col = "blue", lwd = 2)

## Prepare the data for `inla.stack` as well as the projection matrices that
## will be required. Information about vertices will be first, then information
## about observed points.

## Observations
y_pp <- rep(0:1, c(n_vert, n_pp))
## Exposures (expected in each dual mesh polygon)
exposure_pp <- c(wts, rep(0, n_pp))

## Projection matrix to each vertex (Diagonal)
vert_proj <- Diagonal(n_vert, rep(1, n_vert))
## Projection matrix to each observed point
obs_proj <- inla.spde.make.A(mesh, pp_mat)
## Total projection matrix
proj_mat <- rbind(vert_proj, obs_proj)

## Form the inla "stack", which contains data and projection matrices in ways
## that INLA understands.
pp_stack <- inla.stack(data = list(y = y_pp, e = exposure_pp),
                       A = list(1, proj_mat),
                       effects = list(list(b0 = rep(1, n_vert + n_pp)),
                                      list(i = 1:n_vert)),
                       tag = "pp")

###=============================================================================
### Model fitting---------------------------------------------------------------

## ## First fit the model without the areal observations
pp_fit <- inla(y ~ 0 + b0 + f(i, model = spde),
               family = "poisson",
               data = inla.stack.data(pp_stack),
               control.predictor = list(A = inla.stack.A(pp_stack),
                                        compute = TRUE),
               E = inla.stack.data(pp_stack)$e)

cmp <- coordinates ~ -1 + intercept + pp_spat(map = coordinates, model = spde)

lgcp_fit <- lgcp(cmp, pp_df,
                 samplers = pp_dom_poly,
                 domain = list(coordinates = mesh))

pxl <- pixels(mesh, nx = 128, ny = 128, mask = pp_dom_poly)

lgcp_pred <- predict(lgcp_fit,
                     pixels(mesh, mask = pp_dom_poly),
                     ~ exp(intercept + pp_spat))

ggplot() +
  gg(lgcp_pred) +
  gg(pp_dom_poly) +
  ggtitle("LGCP fit to point process", subtitle = "(Intensity surface)") +
  coord_fixed() +
  scale_fill_viridis(option = "A")
