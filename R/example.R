# define an ACV function
acv <- function(theta, tau) {
  A <- abs(theta[1])
  l <- abs(theta[2])
  acov <- A * exp(-(tau / l))
  return(acov)
}

# define parameters of ACV
# theta[1] = mu (mean)
# theta[2] = nu (error scale factor)
# theta[3:p] = parameters of ACV
theta <- c(0.0, 1.0, 1.0, 1.0)

# define vector of times for observations
m <- 1000
t <- seq(-0.5, 10.5, length = m)

# define parameters of ACV
# theta[1] = mu (mean)
# theta[2] = nu (error scale factor)
# theta[3:p] = parameters of ACV
theta <- c(0.0, 1.0, 1.0, 5)

# produce Gaussian vector (simulation)
y <- gp_sim(theta, acv.model = acv, t.star = t)
y <- as.vector(y)

# plot the 'true' curve
dat.true <- data.frame(t = t, y = y)
plot_ts(dat.true, type = "l")

# --------------------------------------------
# examine some data

plot_ts(drw, extend = 1.3)

# fit the model to the data, find Max.Like parameter values
result <- gp_fit(theta, acv.model = acv, dat = drw)

# define vector of times for reconstruction
m <- 1000
t <- seq(-0.5, 10.5, length = m)

# reconstruct process: compute 'conditional' mean and covariance
gp <- gp_conditional(result$par, acv.model = acv, drw, t.star = t)

# plot a 'snake' showing mean +/- std.dev
plot_snake(gp, add = TRUE, col.line = "red", col.fill = "pink", sigma = 1)
