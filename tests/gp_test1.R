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

# --------------------------------------
# now select only a subset of random times to make
# observations

n <- 40
indx <- sample(50:950, size = n)
indx <- sort(indx)

t <- t[indx]
y <- y[indx]

# now add measurement errors
dy <- rep(0.05, n)
y <- y + rnorm(n, mean = 0, sd = dy)

# package into data.frame
dat <- data.frame(t = t, y = y, dy = dy)

# --------------------------------------
# plot the observations
plot_ts(dat, add = TRUE, col = "red", type = "p")

# --------------------------------------
# add some reconstrution realisations
t.star <- seq(-1, 12, length = 400)
sims <- gp_sim(theta, dat, t.star = t.star, N.sim = 10,
               plot = TRUE, acv.model = acv)

# --------------------------------------
# compute likelihood

# fit the model to the data, find Max.Like parameter values
result <- gp_fit(theta, acv.model = acv, dat = drw,
                 method = "Nelder-Mead",
                 maxit = 1000, chatter=1)

# --------------------------------------
# extract the parameters of the ACV model
result$par[2:4] <- abs(result$par[2:4])

# compute 'conditional' mean and covariance
gp <- gp_conditional(theta, acv.model = acv, dat, t.star)

# plot a 'snake' showing mean +/- std.dev
plot_snake(gp, add = TRUE, col.line = "red", col.fill = "pink", sigma = 2)

# plot the observations
points(t, y, pch = 16, cex = 1.5)

# plot error bars
segments(t, y-dy, t, y+dy)

# --------------------------------------
# Now use MCMC to estimate posterior density of parameters

# define the 'priors' for the parameter values
logPrior <- function(theta) {
  mu.d <- dnorm(theta[1], sd = 5, log = TRUE)
  nu.d <- dlnorm(theta[2], meanlog = 0, sdlog = 1, log = TRUE)
  A.d <- dlnorm(theta[3], meanlog = 0, sdlog = 1, log = TRUE)
  l.d <- dlnorm(theta[4], meanlog = 0, sdlog = 1, log = TRUE)
  return(mu.d + nu.d + A.d + l.d)
}

# Use gw.mcmc to generate parameter samples
chain <- tonic::gw_sampler(gp_logPosterior, theta.0 = theta,
                           acv.model = acv, logPrior = logPrior,
                           dat = dat, burn.in = 1e4,
                           nsteps = 20e4,
                           chatter = 1, thin = 10)

# name the parameters
colnames(chain$theta) <- c("mu", "nu", "A", "l")

# plot scatter diagrams
tonic::contour_matrix(chain$theta, prob.levels = c(1,2,3), sigma = TRUE,
                      prob1d = 0)

# plot MCMC diagnostics
tonic::chain_diagnosis(chain)

# present posterior inferences
summary(coda::mcmc(chain$theta))

