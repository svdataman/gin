# define covariance function
acv <- function(theta, tau) {
  A <- abs(theta[1])
  l <- abs(theta[2])
  acov <- A * exp(-(tau / l))
  return(acov)
}

# --------------------------------------
#set.seed(733765)

# define vector of times for observations
n <- 25
t <- runif(n, min = 0, max = 1)
# sort into order
t <- sort(t)

theta <- c(1.0, 0.2)
mu <- array(0, dim = n) # set all means to zero
tau <- outer(t, t, "-") # compute t_j - t_i
tau <- abs(tau) # compute |t_j - t_i|
S <- acv(theta, tau) # acf(tau)
diag(S) <- diag(S) + 0.001

# --------------------------------------
# produce Gaussian vector
y <- mvtnorm::rmvnorm(1, mean = mu, sigma = S, method = "chol")
y <- as.vector(y)
# now add measurement errors
dy <- rep(0.05, n)
y <- y + rnorm(n, mean = 0, sd = dy)
# plot the observations
plot(t, y, bty = "n", xlab = "time", ylab = "y", xlim = c(-0.1, 1.1))
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

# no 'predict' the underlying GP y(t) at times t.star
gp <- predict.gp(theta, dat, t.star)

# plot the 'snake of mean +/- sd
plot.snake(gp, add = TRUE)

# add the data points
points(t, y)

# add the error bars
segments(t, y-dy, t, y+dy)

# --------------------------------------
# add some simulations
# sim.gp(theta, dat, t.star = t.star, N.sim = 5, plot = TRUE)

# --------------------------------------
# compute likelihood
#l <- loglike(theta, dat = dat)

theta <- c(0, 1.0, 1.0, 0.3)
result <- fit.gp(theta, dat = dat, method = "Nelder-Mead", maxit = 1000, chatter=1)

result$par[2:4] <- abs(result$par[2:4]) 

gp <- predict.gp(result$par, dat, t.star)
plot.snake(gp, add = TRUE, col.line = "blue", col.snake = col.b)

