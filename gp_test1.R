# define a auto-covariance function (ACV)
acv <- function(theta, tau) {
  A <- abs(theta[1])
  l <- abs(theta[2])
  acov <- A * exp(-(tau / l))
  return(acov)
}

# --------------------------------------
# set.seed(733765)

# define vector of times for observations
m <- 1000
t <- seq(-0.5, 1.5, length = m)

# define parameters of ACV
theta <- c(1.0, 0.2)
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
dy <- rep(0.1, n)
y <- y + rnorm(n, mean = 0, sd = dy)

# plot the observations
points(t, y, pch = 16, cex = 1.5)

# plot error bars
segments(t, y-dy, t, y+dy)

# --------------------------------------
# now try reconstructing 
col.b <- rgb(173, 216, 230, 100, maxColorValue = 255)

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

theta <- c(0, 1.0, 1.0, 0.3)
result <- fit.gp(theta, dat = dat, 
                 method = "Nelder-Mead", 
                 maxit = 1000, chatter=1)

result$par[2:4] <- abs(result$par[2:4]) 

gp <- predict.gp(result$par, dat, t.star)
plot.snake(gp, add = TRUE, col.line = "blue", 
           col.snake = col.b, sigma = 2)

# plot the observations
points(t, y, pch = 16, cex = 1.5)

# plot error bars
segments(t, y-dy, t, y+dy)
