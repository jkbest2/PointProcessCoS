library(INLA)
library(TMB)
library(spatstat)
library(RandomFields)
library(tidyverse)

set.seed(12345)

### Simulate data --------------------------------------------------------------
xrange <- c(0, 128)
yrange <- c(0, 128)
ppwin <- owin(xrange = xrange, yrange = yrange)

pp1 <- rLGCP(model = "matern", nu = 1,
             mu = -3.5, var = 1, scale = 50,
             win = ppwin, saveLambda = TRUE)
pp_locs <- data.frame(s1 = pp1$x, s2 = pp1$y)
pp_locmat <- as.matrix(pp_locs)
## Matrix::image(attr(pp1, "Lambda"))
## points(pp1)



make_grid_locs <- function(step, xrange, yrange) {
  hastep <- step / 2
  expand.grid(s1 = seq(xrange[1] + hastep, xrange[2] - hastep, step),
              s2 = seq(yrange[1] + hastep, yrange[2] - hastep, step))
}

## Generate an INLA mesh
step <- 10
hastep <- step / 2
grid_locs <- make_grid_locs(step, xrange, yrange)
mesh <- inla.mesh.2d(loc = grid_locs,
                     max.edge = c(4, 10), offset = c(hastep, 20))
                     ## boundary = domain_bound)
spde_fem <- inla.mesh.fem(mesh, order = 2)

fine_step <- 1
finegrid_area <- fine_step ^ 2
finegrid_locs <- make_grid_locs(fine_step, xrange, yrange)
## Observation matrix
A_fine <- inla.spde.make.A(mesh, loc = as.matrix(finegrid_locs))

## Calculate number of points within fine grid cells
pp_counts <- pp_locs %>%
  mutate(s1_cell = s1 %/% fine_step,
         s2_cell = s2 %/% fine_step) %>%
  group_by(s1_cell, s2_cell) %>%
  summarize(count = n()) %>%
  ungroup() %>%
  right_join(mutate(finegrid_locs,
                    s1_cell = s1 %/% fine_step,
                    s2_cell = s2 %/% fine_step),
             by = c("s1_cell", "s2_cell")) %>%
  mutate(# s1_center = s1_cell * w_s1 + w_s1 / 2,
         # s2_center = s2_cell * w_s2 + w_s2 / 2,
    count = coalesce(count, 0L)) %>%
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
data <- list(N_quads = nrow(obs_small),
             N_vertices = mesh$n,
             count_full = pp1$n,
             count_quads = obs_small$count,
             C0 = spde_fem$c0,
             G1 = spde_fem$g1,
             G2 = spde_fem$g2,
             A_full = A_full,
             A_obs = Reduce(cbind, obs_small$rows))
pars <- list(mu = -4,
             spat = rep(0, mesh$n),
             log_tau = log(120),
             log_kappa2 = log(100))

## Compile and load the model
compile("ppcos.cpp")
dyn.load(dynlib("ppcos"))

## Make a function object
obj <- MakeADFun(data, pars, random = "spat", DLL="ppcos")

## Call function minimizer
fit <- optim(obj$par, obj$fn, obj$gr, method = "BFGS")

## Get parameter uncertainties and convergence diagnostics
rep <- obj$report()
sdr <- sdreport(obj)
## sdr
