
# -----------------------------------------------------------
#
# History:
#  16/12/13 - v0.1 - First working version
#  04/07/16 - v0.2 - Major ewrite
#
# Simon Vaughan, University of Leicester
# -----------------------------------------------------------

#' Compute log likelihood function for Gaussian Process model.
#'
#' \code{gp_logLikelihood} returns the log likelihood for a GP model.
#'
#' @param theta     (vector) parameters for covariance function
#'                   the first element is the mean value mu
#' @param acv.model (name) name of the function to compute ACV(tau|theta)
#' @param tau       (matrix) N*N array of lags at which to compute ACF
#' @param dat       (data frame) an N * 3 data frame/array, 3 columns
#                    give the times, measurement and errors of
#                    the n data points.
#' @param PDcheck   (logical) use Matrix::nearPD to coerse the matrix
#                     C to be positive definite
#' @param chatter   (integer) higher values give more run-time feedback
#'
#' @return
#'  scalar value of log[likelihood(theta)]
#'
#' @section Notes:
#'  Compute the log likelihood for Gaussian Process model with parameters theta
#'  given data \eqn{\{t, y, dy\}}
#'  and an (optional) N*N matrix of lags, tau. See algorithm 2.1 of Rasmussen &
#'  Williams (2006). The input data frame 'dat' should contain three columns:
#'  \code{t}, \code{y}, \code{dy}. \code{t[i]} and \code{y[i]} give the times
#'  and the measured values of those
#'  times. \code{dy} gives the 'error' on the measurements \code{y}, assumed to be
#'  independent Gaussian errors wih standard deviation \code{dy}. If \code{dy}
#'  is not present
#'  we assumine \code{dy[i] = 0} for all \code{i}. The columns \code{t},
#'  \code{y}, and \code{dy} are all \code{n}-element vectors.
#'
#'  For multivariate normal distribution the likelihood is
#'
#'  \deqn{L(\theta) = (2\pi)^{-N/2} * det(C)^{-1/2} * exp(-1/2 *
#'  (y-\mu)^T C^{-1} (y-\mu))}
#'
#'  where y is an N-element vector of data and C is an N*N covariance matrix
#'  (positive, symmetric, semi-definite). We compute \eqn{l = log(L)} which can be
#'  written:
#'
#'  \deqn{  l(\theta) = -(n/2) * log(2\pi)
#'        - (1/2) log(det(C))
#'        - (1/2) ((y-mu)^T C^{-1} (y-mu)}
#'
#'  The N*N matrix inverse \eqn{C^{-1}} is slow. Cholesky decomposition allows a
#'  faster calculation of l:
#'
#'  \deqn{  C = LL^T }
#'
#'  and so
#'
#'  \deqn{  \det(C) = \prod L_{ii}^2 }
#'  \deqn{  \log(\det(C)) = 2 \sum \log L_{ii} = 2 \sum diag(L) }
#'
#'  and
#'
#'  \deqn{  C^{-1} = (L L^T)^{-1} = (L^{-1})^T (L^{-1}) }
#'  \deqn{   Q = (y-mu)^T C^{-1} (y-mu) }
#'  \deqn{     = (y')^T (L^{-1})^T (L^{-1}) (y') }
#'  \deqn{     = [(L^{-1}) (y')]^T [(L^{-1}) (y')] }
#'  \deqn{     = z^T z}
#'   where \eqn{z = L^{-1} y'}, and \eqn{y' = Lz}. We can find \eqn{z} using
#'   \eqn{z = solve(L,y')} where \eqn{y' = y - mu}.
#'
#' @seealso \code{\link{gp_logPosterior}}
#'
#' @export
gp_logLikelihood <- function(theta,
                    acv.model = NULL,
                    tau = NULL,
                    dat = NULL,
                    PDcheck = TRUE,
                    chatter = 0) {

  # check arguments
  if (missing(theta)) {stop('** Missing theta input.')}
  if (!all(is.finite(theta))) {stop('** Non-finite values in theta.')}
  if (is.null(dat)) {stop('** Missing dat.')}
  if (!("y" %in% names(dat))) {stop('** Missing dat$y.')}
  if (!("t" %in% names(dat))) {stop('** Missing dat$t.')}
  if (missing(acv.model)) stop('Must specify name of ACV function')
  if (!exists('acv.model')) stop('The specified ACV function does not exist.')

  # length of data vector(s)
  n <- length(dat$y)

  # if there are no errors, and the dat$dy column is
  # missing, make a column of zeroes.
  if (!("dy" %in% colnames(dat))) {
    dat$dy <- array(0, n)
  }
  dy <- dat$dy

  # if n * n array tau is not present then make one
  if (is.null(tau)) {
    tau <- gp_lagMatrix(dat$t)
  }

  # make sure y and tau have matching lengths
  if (ncol(tau) != n) {
    stop('** y and tau do not have matching dimensions.')
  }

  # first, extract the mean (mu) from the parameter vector theta,
  # the extract the error scaling parameter (nu) from the vector theta,
  # Then remove these the theta, so now these only contains parameter for
  # the ACV function
  if (chatter > 1) { cat(theta) }
  mu <- theta[1]
  nu <- abs(theta[2])
  theta <- theta[c(-1, -2)]

  # now subtract the mean from the data: y <- (y - mu)
  y <- dat$y - mu

  # compute the covariance matrix C as C[i,j] = ACV(tau[i,j])
  # using the remaining parameters
  C <- acv.model(theta, tau)

  # check there aren't any non-finite values. If there are then we set the log
  # likelihood to = -Inf.
  if (!all(is.finite(C))) {
    cat('Non-finite values in model covariance matrix.')
    return(-Inf)
  }

  # add the error matrix diag(dy^2)
  diag(C) <- diag(C) + nu*dy*dy

  # enforce for positive-definiteness (PD) and symmetry The covariance matrix
  # should be a PD matrix but numerical errors may mean this is not exactly
  # true. We use the nearPD function from the Matrix package to find the nearest
  # PD matrix and use his.
  if (PDcheck == TRUE) {
    pd <- Matrix::nearPD(C)
    if (pd$converged == FALSE) {
      warning('** nearPD failed to converge.')
    }
    C <- pd$mat
    C <- matrix(C, nrow = n)
    rm(pd)
  }

  # compute the easy (constant) part of the log likelihood.
  l.1 <- -(n/2) * log(2*pi)

  # find the Cholesky decomposition (lower left)
  L <- t( chol(C) )

  # first, compute the log|C| term
  l.2 <- -sum( log( diag(L) ) )

  # then compute the quadratic form -(1/2) * (y-mu)^T C^-1 (y-mu)
  z <- as.vector( solve(L, y) )
  l.3 <- -0.5 * (z %*% z)[1]

  # combine all three terms for give the log[likelihood]
  llike <- l.1 + l.2 + l.3
  if (chatter > 1) { cat(' ', llike, fill = TRUE) }
  return(llike)
}

# -----------------------------------------------------------

#' Compute a posterior density given likelihood and prior.
#'
#' \code{gp_logPosterior} compute posterior density for GP model.
#'
#' Calculate posterior density for a Gaussian Process (GP) model
#' given functions for the AutoCovariance (ACV) of the GP, the
#' log likelihood, and the log prior (both log densities).
#'
#' @param PDcheck (logical) check/correct for positive definiteness?
#' @inheritParams gp_logLikelihood
#'
#' @return
#'   log posterior density (scalar)
#'
#' @seealso \code{\link{gp_logLikelihood}}
#'
#' @export
gp_logPosterior <- function(theta,
                         acv.model = NULL,
                         tau = NULL,
                         dat = NULL,
                         PDcheck = TRUE,
                         chatter = 0,
                         logPrior = NULL) {

  # check arguments
  if (missing(theta)) {stop('** Missing theta input.')}
  if (!all(is.finite(theta))) {stop('** Non-finite values in theta.')}
  if (is.null(dat)) {stop('** Missing dat.')}
  if (!("y" %in% names(dat))) {stop('** Missing dat$y.')}
  if (!("t" %in% names(dat))) {stop('** Missing dat$t.')}
  if (missing(acv.model)) stop('Must specify name of ACV function')
  if (!exists('acv.model')) stop('The specified ACV function does not exist.')

  if (!exists('gp_logLikelihood')) stop('The gp_logLikelihood function does not exist.')

  # compute prior
  if (is.null(logPrior)) {
    lprior <- 0
  } else {
    if (!exists('logPrior')) stop('The logPrior function does not exist.')
    lprior <- logPrior(theta)
  }

  # check if prior = 0; if so, no need to compuite gp_logLikelihood
  if (lprior == -Inf) {
    llike <- 0
  } else {
    llike <- gp_logLikelihood(theta, acv.model, tau = tau,
                   dat = dat, PDcheck = PDcheck,
                   chatter = chatter)
  }

  # posterior ~ likelihood * prior
  lpost <- llike + lprior
  return(lpost)
}

# -----------------------------------------------------------
# Build the matrix of lags tau[i, j] = |t.i - t.j|
#
# (c) Simon Vaughan, University of Leicester
# -----------------------------------------------------------

#' Compute a matrix of differences given two vectors.
#'
#' \code{gp_lagMatrix} returns a matrix of differences between two vectors.
#'
#' Given two vectors - \code{x} (length \code{M}) and \code{y} (length \code{N}) - as
#' input, return the \code{N*M} matrix of differences \code{result[i,j] = x[i] - y[j]}.
#'
#' @param x vector 1
#' @param y vector 2 (default to vector 1 if not specified)
#'
#' @return
#' \code{N*M} array of differences, \code{result[i,j] = x[i] - y[j]}
#'
#' @section Notes:
#' Note that in the special case that \code{x=y} we have a square symmetric
#' matrix: \code{result[i,j] = result[j,i]}. In the even more special case that
#' the two vectors are evenly spaced (\code{x[i] = y[i] = i * delta + const})
#' then we have a circulant matrix; the \code{j}th column \code{result[,j]} is
#' the \code{(j-1)}th cyclic permutation of the first column. This matrix is
#' symmetric, Toeplitz and circulant.
#'
#' @examples
#' result <- gp_lagMatrix(c(1,2,3), c(2,3,4,5,6))
#' print(result)
#'
#' @export
gp_lagMatrix <- function(x, y) {

  # check arguments
  if (missing(x)) stop('Missing input vector.')
  if (missing(y)) y <- x

  # compute x.0, an arbitrary start point for vector x.
  # This improves accuracy if the offset is large
  # compared to the range of values
  x.0 <- min(x, y)
  x <- x - x.0
  y <- y - x.0

  # compute the time lag matrix
  tau <- abs(outer(x, y, "-"))

  # return to calling function
  return(as.matrix(tau))
}

# -----------------------------------------------------------
# optimise the deviance (-2*log[likelihood])

gp_fit <- function(theta.0,
                   acv.model = NULL,
                   dat = NULL,
                   method = "Nelder-Mead",
                   trace = 0,
                   theta.scale = NULL,
                   maxit = 5E3,
                   chatter = 0,
                   PDcheck = FALSE) {

  # -----------------------------------------------------------
  # Inputs:
  #   theta  - vector of (hyper-)parameters for ACV/PSD
  #   acv.model - name of the function to compute ACV(tau|theta)
  #   dat    - 3 column data frame (or list) containing
  #     t    - vector of n observation times
  #     y    - vector of n observations
  #     dy   - vector of n 'errors' on observations
  #   method - choice of method for optim()
  #   theta.scale - vector of rescaling values for the
  #                 parameters these. optim() will fit
  #                 theta.0/theta.scale.
  #   maxit  - max number of interations, passed to optim()
  #   PDcheck   - TRUE/FALSE use Matrix::nearPD to coerse the matrix
  #                C to be positive definite (passed to gp_logLikelihood)
  #   chatter   - higher values give more run-time feedback
  #
  # Value:
  #  a list containing
  #   par             - parameter values (maximum likelihood estimates)
  #   err             - std. dev. of MLEs (based on Hessian matrix)
  #   acv.model       - name of function used to compute ACV
  #   value           - value of gp_logLikelihood at maximum
  #   covergence      - convergence output from optim
  #   nfunction.calls - counts output from optim()
  #
  # Description:
  #  Find the Maximum Likelihood Estimates (MLEs) for the
  #  parameters theta of the ACF/PSD model, given the data
  #  {x, y, dy}.
  #
  # History:
  #  16/12/13 - v0.1 - First working version
  #  05/07/16 - v0.2 - Major re-write
  #
  # Simon Vaughan, University of Leicester
  # -----------------------------------------------------------

  # check arguments
  if (missing(theta.0)) {stop('** Missing theta.0 input.')}
  if (!all(is.finite(theta.0))) {stop('** Non-finite values in theta.')}
  if (is.null(dat)) {stop('** Missing dat.')}
  if (!("y" %in% names(dat))) {stop('** Missing dat$y.')}
  if (!("t" %in% names(dat))) {stop('** Missing dat$t.')}
  if (is.null(theta.scale)) { theta.scale <- rep(1, length(theta.0)) }
  if (missing(acv.model)) stop('Must specify name of ACV function')
  if (!exists('acv.model')) stop('The specified ACV function does not exist.')

  # no. parameters?
  n.parm <- length(theta.0)

  # compute all the time differences: tau[i,j] = |t[j] - t[i]|
  tau <- gp_lagMatrix(dat$t)

  # check the initial position
  post.0 <- gp_logPosterior(theta.0, tau = tau, dat = dat,
                            acv.model = acv.model)
  if (!is.finite(post.0)) {
    stop('Non-finite log posterior for initial theta value.')
  }

  # perform the fitting
  result <- optim(fn = gp_logPosterior,
                  par = theta.0,
                  method = method,
                  tau = tau,
                  dat = dat,
                  PDcheck = PDcheck,
                  acv.model = acv.model,
                  hessian=TRUE,
                  control=list(fnscale = -abs(post.0),
                               trace = trace,
                               parscale = theta.scale,
                               maxit = maxit))

  # test if Hessian is non-singular. If so, estimate errors using covariance
  # matrix. Otherwise, set errors = 0. The covariance matrix = Inv[ -Hessian ]
  # where Hessian_{ij} = dL^2 / dtheta_i dtheta_j
  if (class(try(solve(-result$hessian),silent = T)) == "matrix") {
    covar <- solve(-result$hessian)
    err <- sqrt(diag(covar))
  } else {
    err <- array(NA, dim = n.parm)
  }

  # output the MLEs to screen.
  if (chatter > 0) {
    for (i in 1:n.parm) {
      cat('-- Parameter ', i, ': ', signif(result$par[i], 4), ' +/- ',
          signif(err[i], 3), fill=TRUE, sep='')
    }
  }

  # return the parameter values and their errors
  outp <- list(par = result$par, err = err,
               acv.model = deparse(substitute(acv.model)),
               value = result$value,
               convergence = result$convergence,
               nfunction.calls = result$counts)
  return(outp)
}

# -----------------------------------------------------------
# Predict the mean of the Gaussian Process

gp_conditional <- function(theta,
                       acv.model = NULL,
                       dat = NULL,
                       t.star = NULL,
                       PDcheck = FALSE) {

  # -----------------------------------------------------------
  # gp_conditional
  # Inputs:
  #   theta     - vector of (hyper-)parameters for ACV/PSD
  #   acv.model - name of the function to compute ACV(tau|theta)
  #   dat       - 3 column data frame (or list) containing
  #     t       - vector of n observation times
  #     y       - vector of n observations
  #     dy      - vector of n 'errors' on observations
  #  t.star     - vector of m times at which to predict
  #   PDcheck   - TRUE/FALSE use Matrix::nearPD to coerse the matrix
  #                 C to be positive definite
  #
  # Value:
  #  result - list containing
  #     t   - prediction times
  #     y   - mean of GP model at times t
  #    dy   - sqrt(variance) of GP model at times t
  #   cov   - full m*m covariance matrix
  #
  # Description:
  # Predict the expectation of the GP with covariance matrix specified by
  # (hyper-)parameters 'theta' at times t.star based on observations y.obs and
  # times t.obs. Uses eqn 2.23 of Rasmussen & Williams (2006).
  #
  # Computes E[y(t.star)] given y(t.obs) and theta, and also compute covariances
  # for times t.star C[i,j] = C(t.star[i], t.star[j]). This may be a large
  # matrix. Also return dy = sqrt( diag( cov ) ).
  #
  # Use the following matricies each defined as
  # K[t1[i],t2[j]] = ACF(|t2[j] - t1[i]|)
  # if we have
  #  t.obs - times at observations
  #  t.star  - times at predictions
  #  K     - ACV(t.obs, t.obs)
  #  K.ik  - ACV(t.obs, t.star)
  #  K.ki  - ACV(t.star, t.obs)
  #  K.kk  - ACV(t.star, t.star)
  #
  # History:
  #  16/12/13 - v0.1 - First working version
  #  05/07/16 - v0.2 - Major re-write
  #
  # Simon Vaughan, University of Leicester
  # -----------------------------------------------------------

  # check arguments
  if (missing(theta)) {stop('** Missing theta argument.')}
  if (!all(is.finite(theta))) {stop('** Non-finite values in theta.')}
  if (is.null(dat)) {stop('** Missing dat.')}
  if (!("y" %in% names(dat))) {stop('** Missing dat$y.')}
  if (!("t" %in% names(dat))) {stop('** Missing dat$t.')}
  if (missing(acv.model)) stop('Must specify name of ACV function')
  if (!exists('acv.model')) stop('The specified ACV function does not exist.')

  # if times of predictions are not given, use times of the observations
  if (is.null(t.star)) { t.star <- dat$t }

  # first, extract the mean (mu) from the parameter vector theta,
  # the extract the error scaling parameter (nu) from the vector theta,
  # Then remove these the theta, so now theta only contains parameter for
  # the ACV function
  mu <- theta[1]
  nu <- theta[2]
  theta <- theta[c(-1, -2)]

  # compute model covariance matrix at observed delays
  tau.obs <- gp_lagMatrix(dat$t, dat$t)
  K <- acv.model(theta, tau.obs)

  # compute matrix of covariances between observations and predictions
  tau.ki <- gp_lagMatrix(t.star, dat$t)
  tau.kk <- gp_lagMatrix(t.star, t.star)
  K.ki <- acv.model(theta, tau.ki)
  K.kk <- acv.model(theta, tau.kk)

  # clean up memory
  rm(tau.obs, tau.ki, tau.kk)

  # eqn 2.20 of R&W - add the "error" term to the covariance matrix
  C <- K + nu * dat$dy^2 * diag(1, NCOL(K))
  rm(K)

  # compute the inverse covariance matrix
  C.inv <- solve(C)

  # eqn 2.23 of R&W - predict the mean value at prediction times
  y.star <- as.vector(K.ki %*% C.inv %*% (dat$y - mu)) + mu

  # eqn 2.24 of R&W - predict the covariances at all prediction times
  cov.star <- K.kk - K.ki %*% C.inv %*% t(K.ki)

  # enforce positive-definiteness (PD) and symmetry of the resulting
  # covariance matrix. This should be an m*m PD matrix but numerical errors may
  # mean this is not exactly true. We use the nearPD function from the Matrix
  # package to find the nearest PD matrix and use his.
  if (PDcheck == TRUE) {
    pd <- Matrix::nearPD(cov.star)
    if (pd$converged == FALSE) {
      warning('** nearPD failed to converge.')
    }
    cov.star <- pd$mat
    cov.star <- as.matrix(cov.star)
    rm(pd)
  }

  # Extract the diagonal elements from the covariance matrix, i.e. the variances
  # at times t.star. Store this as a vector. Use the absolute value to avoid any
  # negative values creeping in due to numerical errors.
  dy.star <- sqrt( as.vector( abs( diag(cov.star) ) ) )

  # define the output product
  result <- list(t = t.star, y = y.star, dy = dy.star, cov = cov.star)
  return(result)
}

# -----------------------------------------------------------
# generate a random GP realisation

gp_sim <- function(theta,
                   acv.model = NULL,
                   dat = NULL,
                   t.star = NULL,
                   N.sim = 1,
                   plot = FALSE) {

  # -----------------------------------------------------------
  # gp_sim
  # Inputs:
  #  theta   - vector of (hyper-)parameters for ACV/PSD
  #  acv.model - name of the function to compute ACV(tau|theta)
  #  dat     - 3 column data frame (or list) containing
  #     t    - vector of n observation times
  #     y    - vector of n observations
  #     dy   - vector of n 'errors' on observations
  #  t.star  - vector of m times at which to simulate
  #  N.sim   - how many simulations to produce
  #  plot    - TRUE/FALSE - overlay a plot of the simulations?
  #
  # Value:
  #  y.out   - [n, N.sim] matrix, each column a simulation
  #
  # Description:
  # Simulate random realisations of a Gaussian Process (GP) with covariance
  # defined by parameters theta, at times t.star, given data {t, y, dy}. This is
  # essentially eqn 2.22 of Rasmussen & Williams (2006):
  #
  #   y.sim ~ N(mean = y.star, cov = C)
  #
  # Uses the function rmvnorm from the mvtnorm package to generate m-dimensional
  # Gaussian vectors.
  #
  # If dat is supplied (a data frame/list containing t and y) then the
  # conditional mean values at times t.star, y[t.star] and covariance matrix C
  # for lags tau[i,j] = |t.star[j] - t.star[i]| are generated by the gp_conditional
  # function.
  #
  # If dat is not supplied we sample from the 'prior' GP, i.e. assume mean =
  # theta[1] and covariance matrix given by acv.model and theta parameters.
  #
  # History:
  #  17/12/13 - v0.1 - First working version
  #  05/07/16 - v0.2 - Major re-write
  #
  # Simon Vaughan, University of Leicester
  # -----------------------------------------------------------

  # check arguments
  if (missing(theta)) {stop('** Missing theta argument.')}
  if (!all(is.finite(theta))) {stop('** Non-finite values in theta.')}
  if (is.null(dat)) {
    prior.sim <- TRUE
  } else {
    prior.sim <- FALSE
    if (!("y" %in% names(dat))) {stop('** Missing dat$y.')}
    if (!("t" %in% names(dat))) {stop('** Missing dat$t.')}
    n <- length(dat$y)
  }
  if (missing(acv.model)) stop('Must specify name of ACV function')
  if (!exists('acv.model')) stop('The specified ACV function does not exist.')

  # if times of predictions are not given, use times of the observations
  if (is.null(t.star)) {
    if (is.null(dat)) {
      stop('** Must specify dat and/or t.star.')
    }
    t.star <- dat$t
  }

  # length of data to predict/simulated
  m <- length(t.star)

  if (prior.sim == FALSE) {
  # compute the mean and covariance matrix for the GP at times t.star
    gp.out <- gp_conditional(theta, acv.model, dat, t.star)
  } else {
    tau.kk <- gp_lagMatrix(t.star, t.star)
    gp.out <- list(y = rep(theta[1], m),
                   cov = acv.model(theta[-c(1, 2)], tau.kk))
    rm(tau.kk)
  }

  # compute each time series, an m-vector drawn from the
  # multivariate (m-dimensional) Gaussian distribution.
  y.out <- mvtnorm::rmvnorm(n = N.sim, mean = gp.out$y,
                            sigma = gp.out$cov, method = "chol")
  y.out <- t(y.out)

  # plot all the time series
  if (plot == TRUE) {
    for (i in 1:N.sim) {
      lines(t.star, y.out[,i], col = "grey60")
    }
  }

  # return the results
  return(y.out)
}

