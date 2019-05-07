library(INLA)
library(TMB)
library(spatstat)
library(RandomFields)
library(tidyverse)
 
set.seed(2122)

### Simulate data --------------------------------------------------------------
## These ranges are set as `c(0, 128)` because that is the discretization that
## `rLGCP` uses for the Lambda matrix. This way each cell in the lambda matrix
## is area one and comparison with our results later is slightly more
## convenient.
xrange <- c(0, 128)
yrange <- c(0, 128)
ppwin <- owin(xrange = xrange, yrange = yrange)

## Generate a log-Gaussian Cox process realization. Save the intensity surface.
## This intensity surface is what we are ultimately trying to make inference on.
rho = 50
sigma2 = 1 ^ 2
## Want to make the underlying mean high enough that we have lots of points.
## Otherwise our point process data won't be very informative. For mu = -2, we
## expect to see exp(-2) * 128^2 = 2217 points.
mu = 0
## Set nu = 1; this corresponds to alpha = 2 in the SPDE formulation.
pp1 <- rLGCP(model = "matern", nu = 1,
             mu = mu, var = sigma2, scale = rho,
             win = ppwin, saveLambda = TRUE)
pp_locs <- data.frame(s1 = pp1$x, s2 = pp1$y)

## Quick plot of the point process and underlying intensity surface
Matrix::image(attr(pp1, "Lambda"))
points(pp1)

##--Set up finite element mesh, extract required matrices-----------------------
## Make data frame of centers of grid squares `step` apart. E.g.
## `make_grid_locs(1, c(0, 2), c(0, 2))`
make_grid_locs <- function(step, xrange, yrange) {
  hastep <- step / 2
  expand.grid(s1 = seq(xrange[1] + hastep, xrange[2] - hastep, step),
              s2 = seq(yrange[1] + hastep, yrange[2] - hastep, step))
}

## Generate an INLA mesh. Use `step` here to control the number of nodes in your
## mesh; here I'm trying for around 2500 nodes to keep computation time
## reasonable. This grid is just used to seed the nodes for the INLA mesh.
step <- 5
hastep <- step / 2
grid_locs <- make_grid_locs(step, xrange, yrange)
mesh <- inla.mesh.2d(loc = grid_locs,
                     max.edge = c(4, 10), offset = c(hastep, 20))
## This extracts the C0, G1, and G2 matrices so that they can be passed to TMB.
spde_fem <- inla.mesh.fem(mesh, order = 2)

## Generate a fine grid of points. This grid is used to calculate the point
## process likelihood; the point observations are discretized into these cells.
## Note that we have chosen our step so that the area of each cell is 1, but
## this will be passed to TMB and can be changed.
fine_step <- 1
finegrid_area <- fine_step ^ 2
finegrid_locs <- make_grid_locs(fine_step, xrange, yrange)
## Form the "observation matrix". This linearly interpolates from the mesh nodes
## to the observation locations (in this case the centers of each cell in our
## fine grid). This returns a sparse matrix with a row for each observation and
## a column for each mesh node. Each row will have at most three non-zero
## elements.
A_fine <- inla.spde.make.A(mesh, loc = as.matrix(finegrid_locs))

## Discretize the observed points into counts within each fine grid cells. Cells
## with a zero count also contribute to the likelihood, so these are added as
## well.
pp_counts <- pp_locs %>%
  ## Use integer division to get the "cell coordinate" (number of cells from
  ## bottom left cell)
  mutate(s1_cell = s1 %/% fine_step,
         s2_cell = s2 %/% fine_step) %>%
  ## Sum the number of observations in each cell
  group_by(s1_cell, s2_cell) %>%
  summarize(count = n()) %>%
  ungroup() %>%
  ## Join with the data frame containing the centers of each fine grid cell.
  ## Calculate the "cell coordinate" if this data frame on the fly, and join by
  ## them so that the cell coordinates can be associated with each cell.
  right_join(mutate(finegrid_locs,
                    s1_cell = s1 %/% fine_step,
                    s2_cell = s2 %/% fine_step),
             by = c("s1_cell", "s2_cell")) %>%
  ## Missing values are imputed during the `join` operation anywhere there is no
  ## count. In this case we know the count must be zero, so use `coalesce` to
  ## replace `NA`s with `0`.
  mutate(count = coalesce(count, 0L)) %>%
  ## Finally, we only need the coordinates of the centers of each cell and its
  ## observed count, so only return these.
  select(s1_center = s1,
         s2_center = s2,
         count = count)

## ggplot(pp_counts, aes(x = s1_center, y = s2_center, fill = factor(count))) +
##   geom_raster() + coord_equal() +
##   geom_point(mapping = aes(x = s1, y = s2), fill = "black", data = pp_locs)

## filter_locs <- function(lim_list, locs) {
##   (locs$s1 >= lim_list$min_s1) & (locs$s1 <= lim_list$max_s1) &
##   (locs$s2 >= lim_list$min_s2) & (locs$s2 <= lim_list$max_s2)
## }
## count_obs <- function(lim_list, pp) {
##   sum((pp1$x >= lim_list$min_s1) & (pp1$x <= lim_list$max_s1) &
##       (pp1$y >= lim_list$min_s2) & (pp1$y <= lim_list$max_s2))
## }

## areal_divs <- 4L
## areal_wd <- 100 / areal_divs
## n_small <- 5L
## s_small <- sample(0L:(areal_divs ^ 2), n_small, replace = FALSE)
## obs_small <- data_frame(id = s_small,
##                         n_s1 = s_small %/% areal_divs,
##                         n_s2 = s_small %% areal_divs,
##                         obs_lims = map2(n_s1, n_s2,
##                                         ~ list(min_s1 = areal_wd * .x,
##                                                max_s1 = areal_wd * .x + areal_wd,
##                                                min_s2 = areal_wd * .y,
##                                                max_s2 = areal_wd * .y + areal_wd)),
##                         rows_lgl = map(obs_lims, filter_locs,
##                                        locs = loc_fine_df),
##                         rows = map(rows_lgl, which),
##                         A = map(rows, ~ A_full[.x, ]),
##                         count = map_int(obs_lims, count_obs))

## Data and parameters
## data <- list(N_quads = nrow(obs_small),
##              N_vertices = mesh$n,
##              count_full = pp1$n,
##              count_quads = obs_small$count,
##              C0 = spde_fem$c0,
##              G1 = spde_fem$g1,
##              G2 = spde_fem$g2,
##              A_full = A_full,
##              A_obs = Reduce(cbind, obs_small$rows))
## pars <- list(mu = -4,
##              spat = rep(0, mesh$n),
##              log_tau = log(120),
##              log_kappa2 = log(100))

## Form the `list` of data to be passed to TMB.
data <- list(N_vertices = mesh$n,
             N_cells = nrow(pp_counts),
             pp_counts = pp_counts$count,
             pp_cell_area = finegrid_area,
             C0 = spde_fem$c0,
             G1 = spde_fem$g1,
             G2 = spde_fem$g2,
             A = A_fine)
## The SPDE approach uses a different parameterization of range and marginal
## variance parameters. Starting with the generative parameters, transform them
## and pass to TMB as starting values (or fix them) to make fitting easier.
kappa2 <- 8 / rho ^ 2
tau <- sqrt(4 * pi * kappa2 * sigma2)
pars <- list(mu = log(pp1$n / area(ppwin)),
             spat = rep(0, mesh$n),
             log_tau = log(tau),
             log_kappa2 = log(kappa2))
map <- list(log_tau = factor(NA))
            ## log_kappa2 = factor(NA))

## Compile and load the model
compile("ppcos.cpp")
dyn.load(dynlib("ppcos"))

## Make a function object
obj <- MakeADFun(data, pars,
                 map = map, random = "spat",
                 DLL = "ppcos")

## Call function minimizer
fit <- nlminb(obj$par, obj$fn, obj$gr)
## fit <- optim(obj$par, obj$fn, obj$gr, method = "BFGS")
h <- optimHess(fit$par, obj$fn, obj$gr)

## Get parameter uncertainties and convergence diagnostics
rep <- obj$report()
sdr <- sdreport(obj)

## Plot the fitted intensity surface and the generative intensity surface
fit_lambda <- matrix(rep$lambda, )
