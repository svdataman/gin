# --------------------------------------
# test using 5548 data
# tests the function ts_load, plot_ts, plot_snake, gp_conditional

# load the data
filename <- system.file("extdata", "ngc5548_W2.csv", package = "gin")
d <- gin::ts_load(filename, header = TRUE, csv = TRUE, cols=c(2, 4, 5))

# plot the data
gin::plot_ts(d)

# given some ACF parameters... (these are just guesses!)
theta <- c(0.0, 1.0, 3.0, 30)

# define an ACV function
acv <- function(theta, tau) {
  A <- abs(theta[1])
  l <- abs(theta[2])
  acov <- A * exp(-(tau / l^2))
  return(acov)
}

# reconstruct the data
gp <- gin::gp_conditional(theta, acv.model = acv, d, t.star = 56700:56835)

# plot the reconstruction
gin::plot_snake(gp, add = TRUE, col.line = 3, sigma = c(1, 2))

# replot data
gin::plot_ts(d, add = TRUE)

# --------------------------------------

