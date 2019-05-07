library(inlabru)

matern <- inla.spde2.pcmatern(mesh,
                              prior.sigma = c(1, 0.5),
                              prior.range = c(50, 0.5))

cmp <- coordinates ~ mySmooth(map = coordinates, model = matern) + Intercept

dat <- cbind(pp1$x, pp1$y)

fit <- lgcp(cmp, dat)
