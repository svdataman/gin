# define a auto-covariance function (ACV)
acv <- function(theta, tau) {
  A <- abs(theta[1])
  l <- abs(theta[2])
  acov <- A * exp(-(tau / l))
  return(acov)
}

acv2 <- function(theta, tau) {
  A <- abs(theta[1])
  l <- abs(theta[2])
  acov <- A * exp(-(tau / l)^2)
  return(acov)
}

# --------------------------------------
# set.seed(733765)

# define vector of times for observations
m <- 1000
t <- seq(-0.5, 1.5, length = m)

# define parameters of ACV
theta <- c(1.0, 0.1)
mu <- array(0, dim = m) # set all means to zero
tau <- matrix.tau(t) # compute |t_j - t_i|
S <- acv(theta, tau) # acf(tau)
diag(S) <- diag(S) + 0.001

# produce Gaussian vector
y <- mvtnorm::rmvnorm(1, mean = mu, sigma = S, method = "chol")
y <- as.vector(y)

# plot the 'true' curve
plot(t, y, bty = "n", xlab = "time", ylab = "y", 
     type = "l", xlim = c(-0.1, 1.1))

# --------------------------------------
# now select only a subset of random times to make
# observations

n <- 30
indx <- sample(251:750, size = n)
indx <- sort(indx)

t <- t[indx]
y <- y[indx]

# now add measurement errors
dy <- rep(0.05, n)
y <- y + rnorm(n, mean = 0, sd = dy)

# plot the observations
points(t, y, pch = 16, cex = 1.5)

# plot error bars
segments(t, y-dy, t, y+dy)

# --------------------------------------
# now try reconstructing 

dat <- data.frame(t = t, y = y, dy = dy)
t.star <- seq(-1, 2.2, length = 400)

# theta[1] = mu (mean)
# theta[2] = nu (error scale factor)
# theta[3:p] = parameters of ACV
theta <- c(0.0, 1.0, 1.0, 0.2)

# now 'predict' the underlying GP y(t) at times t.star
# gp <- predict.gp(theta, dat, t.star)

# plot the 'snake of mean +/- sd
# plot.snake(gp, add = TRUE)

# --------------------------------------
# add some simulations
# sim.gp(theta, dat, t.star = t.star, N.sim = 5, plot = TRUE)

# --------------------------------------
# compute likelihood
#l <- loglike(theta, dat = dat)

# choose some intitial parameter values
theta <- c(0, 1.0, 1.0, 0.3)

# fit the model to the data, find Max.Like parameter values
result <- gp.fit(theta, acv.model = acv, dat = dat, 
                 method = "Nelder-Mead", 
                 maxit = 1000, chatter=1)

# --------------------------------------
# extract the parameters of the ACV model
result$par[2:4] <- abs(result$par[2:4]) 

# compute 'conditional' mean and covariance
gp <- gp.predict(result$par, acv.model = acv, dat, t.star)

# plot a 'snake' showing mean +/- std.dev
plot.snake(gp, add = TRUE, col.line = "red", 
           sigma = 2, col.fill = "pink")

# plot the observations
points(t, y, pch = 16, cex = 1.5)

# plot error bars
segments(t, y-dy, t, y+dy)

# --------------------------------------
# Now use MCMC to estimate posterior density of parameters
stop()

source('~/R/tonic/gwmcmc.R')
source('~/R/tonic/plot_contour.R')
source('~/R/tonic/diagnostics.R')

# define the 'priors' for the parameter values
LogPrior <- function(theta) {
  mu.d <- dnorm(theta[1], sd = 5, log = TRUE)
  nu.d <- dlnorm(theta[2], meanlog = 0, sdlog = 1, log = TRUE)
  A.d <- dlnorm(theta[3], meanlog = 0, sdlog = 1, log = TRUE)
  l.d <- dlnorm(theta[4], meanlog = 0, sdlog = 1, log = TRUE)
  return(mu.d <- nu.d + A.d + l.d)
}

# Use gw.mcmc to generate parameter samples
chain <- gw.mcmc(LogPosterior, c(0.9, 1.9, 0.5, 0.1), acv.model = acv,
    dat = dat, burn.in = 1e4, nwalkers = 20, nsteps = 1e4, chatter = 1)

# name the parameters
colnames(chain$theta) <- c("mu", "nu", "A", "l")

# plot scatter diagrams
cont.pairs(chain$theta)

# plot MCMC diagnostics
mcmc.diag.plot(chain)

