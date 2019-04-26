
# ------------------------------------------------
# illusrate a Gaussian Process
# Begin with the prior, no data at all.
# Then plot the process conditional on 1 data point.
# Then include more data.

# define an ACV function
acv <- function(theta, tau) {
  A <- abs(theta[1])
  l <- abs(theta[2])
  acov <- A * exp(-(tau / l^2))
  return(acov)
}

# define vector of times for reconstruction
t <- seq(0, 100, by = 1)
m <- length(t)

# define ACF parameters
theta <- c(0.0, 1.0, 9.0, 5)

# make some data
t.obs <- c(50, 20, 30, 35)
y.obs <- gin::gp_sim(theta, acv.model = acv, t.star = t.obs)
y.obs <- as.vector(y.obs)
dy <- rep(0.2, length(t.obs))
y.obs <- y.obs + rnorm(length(t.obs), sd = dy)
obs <- data.frame(t = t.obs, y = y.obs, dy = dy)

# produce Gaussian vector - first with no data, zero mean
y <- rep(0, m)
dy2 <- rep(acv(theta[3:4], 0), m)
dy <- sqrt(dy2)
dat <- data.frame(t = t, y = y, dy = dy)

plot(0, 0, type = "n", ylim = c(-8, 8), xlim = range(t), bty = "n",
     xlab = "time", ylab = "y(t)")
plot_snake(dat, add = TRUE, col.line = 3, sigma = c(1,2))

# add one point,
mask <- c(1)
gp <- gp_conditional(theta, acv.model = acv, obs[mask,], t.star = dat$t)
plot(0, 0, type = "n", ylim = c(-8, 8), xlim = range(t), bty = "n",
     xlab = "time", ylab = "y(t)", axes = FALSE, cex.lab = 2)
grid()
axis(1, lwd=0, cex.axis=1.5)
plot_snake(gp, add = TRUE, col.line = 3, sigma = c(1,2))
points(obs$t[mask], obs$y[mask], pch = 16, cex = 1.2)
segments(obs$t[mask], obs$y[mask]-obs$dy[mask],
         obs$t[mask], obs$y[mask]+obs$dy[mask])

# add two points,
mask <- c(1,2,3,4)
gp <- gp_conditional(theta, acv.model = acv, obs[mask,], t.star = dat$t)
plot(0, 0, type = "n", ylim = c(-8, 8), xlim = range(t), bty = "n",
     xlab = "time", ylab = "y(t)", axes = FALSE, cex.lab = 2)
grid()
axis(1, lwd=0, cex.axis=1.5)
plot_snake(gp, add = TRUE, col.line = 3, sigma = c(1,2))
points(obs$t[mask], obs$y[mask], pch = 16, cex = 1.2)
segments(obs$t[mask], obs$y[mask]-obs$dy[mask],
         obs$t[mask], obs$y[mask]+obs$dy[mask])
