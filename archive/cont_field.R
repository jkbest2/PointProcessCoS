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
## Matrix::image(log(attr(pp1, "Lambda")))
## points(pp1)

field <- log(attr(pp1, "Lambda"))

n_obs <- 1000L
obs_loc_x <- runif(n_obs, 0, 128)
obs_loc_y <- runif(n_obs, 0, 128)

obs_cell <- lapply(1:n_obs, function(i) list(ceiling(obs_loc_x[i]),
                                             ceiling(obs_loc_y[i])))

obs <- sapply(obs_cell,
              function(obs_loc) field[obs_loc[[1]], obs_loc[[2]]] +
                                rnorm(1, 0, 0.1))

obs_df <- data.frame(s1 = obs_loc_x,
                     s2 = obs_loc_y,
                     val = obs)
coordinates(obs_df)  <- ~ s1 + s2

##--Set up SPDE mesh etc--------------------------------------------------------
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

## Calculate `A` matrix to project from mesh nodes to observation locations
A <- inla.spde.make.A(mesh, loc = cbind(obs_loc_x, obs_loc_y))

##--Fit the model in TMB--------------------------------------------------------
compile("cont_field.cpp")
dyn.load(dynlib("cont_field"))

data <- list(N_vertices = mesh$n,
             N_obs = n_obs,
             obs_value = obs,
             C0 = spde_fem$c0,
             G1 = spde_fem$g1,
             G2 = spde_fem$g2,
             A = A)
kappa2 <- 8 / rho ^ 2
tau <- sqrt(4 * pi * kappa2 * sigma2)
parameters <- list(mu = 0,
                   sigma = 0.1,
                   spat = rep(0, mesh$n),
                   log_tau = log(tau),
                   log_kappa2 = log(kappa2))
map <- list(log_tau = factor(NA),
            log_kappa2 = factor(NA),
            sigma = factor(NA))

obj <- MakeADFun(data = data,
                 parameters = parameters,
                 map = map,
                 random = "spat",
                 DLL = "cont_field")

opt <- optim(obj$par, obj$fn, obj$gr)

rep <- obj$report()
sdr <- sdreport(obj)

full_grid <- make_grid_locs(1, c(0, 128), c(0, 128))
A_full <- inla.spde.make.A(mesh, loc = as.matrix(full_grid))
spat <- sdr$par.random
spat_full <- A_full %*% spat

fit_field <- sdr$par.fixed[1] + spat_full@x

## par(mfrow = c(1, 2))
## Matrix::image(rev(1:128), 1:128, matrix(spat_full, nrow = 128))
## image(log(attr(pp1, "Lambda")))

df <- data_frame(s1 = rep(full_grid[, 1], 2),
                 s2 = rep(full_grid[, 2], 2),
                 mod = rep(c("fit", "orig"), each = nrow(full_grid)),
                 val = c(fit_field, c(field$v)))

library(viridis)
ggplot(df, aes(x = s1, y = s2, fill = val)) +
  geom_raster() + scale_fill_viridis() +
  facet_wrap(~ mod)

library(inlabru)

pred_df <- data_frame()
dat_df <- bind_cols()

fit <- bru(val ~ Intercept +
             mySPDE(map = coordinates, model = inla.spde2.matern(mesh, alpha = 2)),
           data = obs_df,
           family = "gaussian")

coordinates(pp_locs)  <- ~ s1 + s2
matern <- inla.spde2.pcmatern(mesh,
                              prior.sigma = c(1, 0.5),
                              prior.range = c(50, 0.5))
cmp <- coordinates ~ my_spde(map = coordinates, model = matern) + Intercept

fit <- lgcp(cmp, pp_locs,
            options = bru.options(control.inla = control.inla(int.strategy = "eb")))
