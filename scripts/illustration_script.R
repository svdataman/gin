# define an ACV function
acv <- function(theta, tau) {
  A <- abs(theta[1])
  l <- abs(theta[2])
  acov <- A * exp(-(tau / l^2))
  return(acov)
}

# --------------------------------------

theta <- c(0.0, 1.0, 1.0, 3)
m <- 500
t <- seq(-0.5, 100.5, length = m)
y <- gp_sim(theta, acv.model = acv, t.star = t)
y <- as.vector(y)
dat <- data.frame(cbind(t, y))

# plot the 'true' curve
xlim <- range(dat$t)
ylim <- range(dat$y)
par(mar=c(6, 6, 2, 2))


mask <- seq(1, 500, by = 1)
plot(dat$t[mask], dat$y[mask], col="pink3", type = "l", pch=1, bty = "n",
     axes=FALSE, lwd = 3, cex = 2,
     xlab="time", ylab="y(t)", xlim = xlim, ylim = ylim, cex.lab = 2)
grid()
axis(1, lwd=0, cex.axis=1.5)

mask <- c(100, 250, 300, 350, 400)
#points(dat$t[mask], dat$y[mask], cex = 2, pch = 1, lwd = 3)

y <- gp_sim(theta, acv.model = acv, t.star = t)
y <- as.vector(y)
dat <- data.frame(t = t, y = y)
lines(dat$t, dat$y, lwd = 3, col = "red3")
#points(dat$t[mask], dat$y[mask], cex = 2, pch = 1, lwd = 3)

y <- gp_sim(theta, acv.model = acv, t.star = t)
y <- as.vector(y)
dat <- data.frame(t = t, y = y)
lines(dat$t, dat$y, lwd = 3, col = "red4")
#points(dat$t[mask], dat$y[mask], cex = 2, pch = 1, lwd = 3)

#segments(dat$t[mask], ylim[1], dat$t[mask], ylim[1]+0.2, lwd=3)

# --------------------------------------

# define parameters of ACV
# theta[1] = mu (mean)
# theta[2] = nu (error scale factor)
# theta[3:p] = parameters of ACV
theta <- c(0.0, 1.0, 1.0, 3)

# define vector of times for reconstruction
m <- 1000
t <- seq(-0.5, 100.5, length = m)

# produce Gaussian vector (simulation)
y <- gp_sim(theta, acv.model = acv, t.star = t)
y <- as.vector(y)
dat <- data.frame(t = t, y = y)

# plot the 'true' curve
xlim <- range(dat$t)
ylim <- range(dat$y)
par(mar=c(6,6,2,2))
plot(dat$t, dat$y, col="pink3", type = "l", bty = "n", axes=FALSE, lwd = 3,
     xlab="time", ylab="y(t)", xlim = xlim, ylim = ylim, cex.lab = 2)
grid()
axis(1, lwd=0, cex.axis=1.5)

# ------------------------------------------------
# -------------------------------------------------


# now observe f(t) only at n random times
n <- 25
bookend <- round(m/4)
midrange <- (bookend):(m-bookend)
indx <- sample(midrange, size = n)
indx <- sort(indx)
t <- dat$t[indx]
y <- dat$y[indx]

# now add measurement errors
dy <- rep(0.1, n) * sample(c(1.0,0.5,3), size = n, replace = TRUE,
                           prob = c(0.8,0.10,0.1))
epsilon <- rnorm(n, mean = 0, sd = dy)
y <- y + epsilon

plot(0, 0, type = "n", bty = "n", axes=FALSE, lwd = 2,
     xlab="time", ylab="y(t)", xlim = xlim, ylim = ylim, cex.lab = 2)
grid()
axis(1, lwd=0, cex.axis=1.5)
lines(dat$t, dat$y, col = "pink3", lwd = 3)

# plot the observations
points(t, y, pch = 16, cex = 2)

# plot error bars
segments(t, y-dy, t, y+dy, lwd=3)

obs <- list(t=t, y=y, dy=dy)

# ------------------------------------------------
# --------------------------------------------
# plot snake
plot(0, 0, type = "n", bty = "n", axes=FALSE, lwd = 2,
     xlab="time", ylab="y(t)", xlim = xlim, ylim = ylim, cex.lab = 2)
grid()
axis(1, lwd=0, cex.axis=1.5)
# reconstruct process: compute 'conditional' mean and covariance
gp <- gp_conditional(theta, acv.model = acv, obs, t.star = dat$t)

# plot a 'snake' showing mean +/- std.dev
plot_snake(gp, add = TRUE, col.line = 3, sigma = c(1, 2))

# plot the observations
points(t, y, pch = 16, cex = 2)

# plot error bars
segments(t, y-dy, t, y+dy, lwd=3)

# ------------------------------------------------
# --------------------------------------------
# plot snake with true data
plot(0, 0, type = "n", bty = "n", axes=FALSE, lwd = 2,
     xlab="time", ylab="y(t)", xlim = xlim, ylim = ylim, cex.lab = 2)
grid()
axis(1, lwd=0, cex.axis=1.5)
# plot a 'snake' showing mean +/- std.dev
plot_snake(gp, add = TRUE, col.line = 3, sigma = c(1, 2))

# add some constrained realisations

y.sim <- gp_sim(theta, dat = obs, acv.model = acv, t.star = dat$t,
                N.sim = 5, plot = FALSE)
#for (i in 1:5) lines(dat$t, y.sim[, i], col = i)
lines(dat$t, dat$y, col = "pink3", lwd = 3)

# plot the observations
points(t, y, pch = 16, cex = 2)

# plot error bars
segments(t, y-dy, t, y+dy, lwd=3)


# ------------------------------------------------
# --------------------------------------------
# zoom - plot snake with true data
plot(0, 0, type = "n", bty = "n", axes=FALSE, lwd = 2,
     xlab="time", ylab="y(t)", xlim = c(60, 80), ylim = ylim, cex.lab = 2)
grid()
axis(1, lwd=0, cex.axis=1.5)
# plot a 'snake' showing mean +/- std.dev
plot_snake(gp, add = TRUE, col.line = 3, sigma = c(1, 2))


# add some constrained realisations

y.sim <- gp_sim(theta, dat = obs, acv.model = acv, t.star = dat$t,
                N.sim = 5, plot = FALSE)
#for (i in 1:5) lines(dat$t, y.sim[, i], col = i)
lines(dat$t, dat$y, col = "pink3", lwd =3)

# plot the observations
points(t, y, pch = 16, cex = 2)

# plot error bars
segments(t, y-dy, t, y+dy, lwd=3)


# ------------------------------------------------
# --------------------------------------------
plot(0, 0, type = "n", bty = "n", axes=FALSE, lwd = 2,
     xlab="time", ylab="y(t)", xlim = xlim, ylim = ylim, cex.lab = 2)
grid()
axis(1, lwd=0, cex.axis=1.5)
# plot a 'snake' showing mean +/- std.dev
plot_snake(gp, add = TRUE, col.line = 3, sigma = c(1, 2))

# add some constrained realisations
t.star <- seq(-0.5, 100.5, length = 500)

y.sim <- gp_sim(theta, dat = obs, acv.model = acv, t.star = t.star,
                N.sim = 5, plot = FALSE)
col <- c("red", "blue", "black", "brown")
col <- rgb(t(col2rgb(col)), alpha = 30,
                maxColorValue = 255)
for (i in 1:4) lines(t.star, y.sim[, i], col = col[i], lwd=2)

# plot the observations
points(t, y, pch = 16, cex = 2)

# plot error bars
segments(t, y-dy, t, y+dy, lwd=3)

# ------------------------------------------------


# ------------------------------------------------
# ------------------------------------------------
# ------------------------------------------------
# ------------------------------------------------
# ------------------------------------------------
