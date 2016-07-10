
# -----------------------------------------------------------

loglike <- function(theta, 
                    acv.model,
                    tau = NULL, 
                    dat = NULL, 
                    PDcheck = TRUE,
                    chatter = 0) {
  
  # -----------------------------------------------------------
  # loglike
  # Inputs: 
  #   theta     - vector of parameters for covariance function
  #                the first element is the mean value mu
  #   acv.model - name of the function to compute ACV(tau|theta)
  #   tau       - N*N array of lags at which to compute ACF
  #   dat       - an n * 3 data frame/array, 3 columns
  #                give the times, measurement and errors of
  #                the n data points.
  #   PDcheck   - TRUE/FALSE use Matrix::nearPD to coerse the matrix
  #                C to be positive definite
  #   chatter   - higher values give more run-time feedback
  #
  # Value:
  #  logl       - value of log[likelihood(theta)]
  #
  # Description:
  #  Compute the log-likelihood for model parameters theta given data {t, y, dy}
  #  and an (optional) n*N matrix of lags, tau. See algorithm 2.1 of Rasmussen &
  #  Williams (2006). The input data frame 'dat' should contain three columns: 
  #  t, y, dy. t[i] and y[i] give the times and the measured values of those 
  #  times. dy gives the 'error' on the measurements y, assumed to be 
  #  independent Gaussian errors wih standard deviation dy. If dy is not present
  #  we assumine dy[i] = 0 for all i. The columns t, y, and dy are all n-element
  #  vectors.
  #  
  #  For multivariate normal distribution the likelihood is
  #  
  #  L(theta) = (2*pi)^(-n/2) * det(C)^(-1/2) * exp(-1/2 *
  #  (y-mu)^T.C^(-1).(y-mu))
  #  
  #  where y is an n-element vector of data and C is an n*n covariance matrix
  #  (positive, symmetric, semi-definite). We compute l = log(L) which can be
  #  written
  #  
  #    l = -(n/2) * log(2*pi) 
  #        - (1/2)*log(det(C)) 
  #        - (1/2) * (y-mu)^T.C^(-1).(y-mu)
  #  
  #  The n*n matrix inverse C^(-1) is slow. Cholesky decomposition allows a
  #  faster calculation of l:
  #  
  #    C = L.L^T 
  #    det(C) = prod( L[i,i]^2 ) 
  #    log(det(C)) = 2 * sum( log(L[i,i])) = 2 * sum(diag(L))
  #
  #  and 
  #
  #    C^(-1) = (L L^T)^(-1) = (L^-1)^T (L^-1)
  #    Q = (y-mu)^T.C^(-1).(y-mu)
  #      = (y')^T.(L^-1)^T.(L^-1).(y')
  #      = [(L^-1).(y')]^T.[(L^-1).(y')]
  #      = z^T.z
  #   where z = L^(-1).y' --> y' = L.z --> z = solve(L,y')
  #   and y' = y - mu
  # 
  # History:
  #  16/12/13 - v0.1 - First working version
  #  04/07/16 - v0.2 - Major ewrite
  #
  # Simon Vaughan, University of Leicester
  # -----------------------------------------------------------

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
  
  # if n * n array tau is not present then make one
  if (is.null(tau)) {
    tau <- matrix.tau(dat$t)
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
  loglike <- l.1 + l.2 + l.3
  if (chatter > 1) { cat(' ', loglike, fill = TRUE) }
  return(loglike) 
}

# -----------------------------------------------------------
# Build the matrix of lags tau[i,j] = |t.i - t.j|

matrix.tau <- function(t.i, t.j) {
  
  # -----------------------------------------------------------
  # matrix.tau
  # Inputs: 
  #   t.i   - times t.1 for vector y[t.1[i]]
  #   t.j   - times t.2 for vector y[t.2[j]]
  #
  # Value:
  #  tau.ij - matrix of time lags 
  #
  # Description:
  #  Define the matrix of time lags
  #   tau[i,j] = |t.1[i] - t.2[j]|
  #
  # Note that in the special case that t.1=t.2 we have
  # a square symmetric matrix: tau[i,j] = tau[j,i]. 
  # In the even more special case that
  # the two series are identically and evenly sampled
  # (t.1[i] = t.2[i] = i * dt + t0) then we have a 
  # circulant matrix; the jth column tau[,j] is the (j-1)th 
  # cyclic permutation of the first column. This matrix is
  # symmetric, Toeplitz and circulant. 
  #
  # History:
  #  16/12/13 - First working version
  #
  # Simon Vaughan, University of Leicester
  # -----------------------------------------------------------
  
  # check arguments
  if (missing(t.i)) { stop('** Missing t.i array.')}
  if (missing(t.j)) { t.j <- t.i }

  # subtract off t.0, the arbitrary start time 
  # to make the |t.i - t.j| calculation easier.
  t.0 <- t.i[1]
  t.i <- t.i - t.0
  t.j <- t.j - t.0

  # compute the time lag matrix
  tau <- abs( outer(t.i, t.j, "-") )

  # return to calling function
  return(tau)  
}

# -----------------------------------------------------------
# optimise the deviance (-2*log[likelihood])

fit.gp <- function(theta.0, 
                   acv.model,
                   dat = NULL, 
                   method="Nelder-Mead", 
                   trace=0, 
                   theta.scale=NULL,
                   maxit = 5E3,
                   chatter = 0,
                   PDcheck = FALSE) {
  
  # -----------------------------------------------------------
  # fig.gp
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
  #                C to be positive definite (passed to loglike)
  #   chatter   - higher values give more run-time feedback
  #
  # Value:
  #  outp - list of parameters (par) and their errors (err)
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
  
  # compute all the time differences 
  #  tau[i,j] = |t[j] - t[i]|
  tau <- matrix.tau(dat$t)
  
  # check the initial position
  loglike.0 <- loglike(theta.0, tau = tau, dat = dat, 
                       acv.model = acv.model)
  if (!is.finite(loglike.0)) {
    stop('Non-finite log likelihood for initial theta value.')
  }
  
  # perform the fitting  
  result <- optim(fn = loglike, 
                  par = theta.0, 
                  method = method,
                  tau = tau, 
                  dat = dat,
                  PDcheck = PDcheck,
                  acv.model = acv.model,
                  hessian=TRUE,
                  control=list(fnscale = -abs(loglike.0), 
                               trace = trace, 
                               parscale = theta.scale,
                               maxit = maxit))
  
  if (chatter > 1) { print( (result) ) }
  
  # test if Hessian is non-singular. If so, estimate errors using
  # covariance matrix. Otherwise, set errors = 0.
  # The covariance matrix = Inv[ -Hessian ] where 
  #   Hessian_{ij} = dL^2 / dtheta_i dtheta_j
  if (class(try(solve(-result$hessian),silent = T)) == "matrix") {
    covar <- solve(-result$hessian)
    err <- sqrt(diag(covar))
  } else {
    err <- array(NA, dim = n.parm)
  }

  if (chatter > 0) {  
    for (i in 1:n.parm) {
      cat('-- Parameter ', i, ': ', signif(result$par[i], 4), ' +/- ', 
          signif(err[i], 3), fill=TRUE, sep='')
    }
  }
  
  # return the parameter values and their errors
  outp <- list(par = result$par, err = err, 
               acv.model = deparse(substitute(acv.model)))
  return(outp)  
}

# -----------------------------------------------------------
# Predict the mean of the Gaussian Process 

predict.gp <- function(theta, 
                       acv.model,
                       dat = NULL, 
                       t.star = NULL,
                       PDcheck = FALSE) {
  
  # -----------------------------------------------------------
  # predict.gp
  # Inputs: 
  #   theta - vector of (hyper-)parameters for ACV/PSD
  #   acv.model - name of the function to compute ACV(tau|theta)
  #   dat    - 3 column data frame (or list) containing
  #     t    - vector of n observation times
  #     y    - vector of n observations
  #     dy   - vector of n 'errors' on observations
  #  t.star  - vector of m times at which to predict 
  #   PDcheck - TRUE/FALSE use Matrix::nearPD to coerse the matrix
  #                C to be positive definite 
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
  # matrix, so we return by default only the leading diagnoal, i.e. the 
  # variances at times t.star.
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
  # Then remove these the theta, so now these only contains parameter for
  # the ACV function
  # extract the mean - mu - from the theta vector
  mu <- theta[1]
  nu <- theta[2]
  theta <- theta[c(-1, -2)]

  # compute model covariance matrix at observed delays
  tau.obs <- matrix.tau(dat$t, dat$t)
  K <- acv.model(theta, tau.obs)
  
  # compute matrix of covariances between observations and predictions
  tau.ki <- matrix.tau(t.star, dat$t)
  tau.kk <- matrix.tau(t.star, t.star)
  K.ki <- acv.model(theta, tau.ki)
  K.kk <- acv.model(theta, tau.kk)
  
  # clean up memory
  rm(tau.obs, tau.ki, tau.kk)

  # eqn 2.20 of R&W - add the "error" term to the covariance matrix
  C <- K + nu * dat$dy^2 * diag(1, NCOL(K))
  
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

sim.gp <- function(theta, 
                   acv.model,
                   dat = NULL, 
                   t.star = NULL, 
                   N.sim = 1, 
                   plot = FALSE) {

  # -----------------------------------------------------------
  # sim.gp
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
  # for lags tau[i,j] = |t.star[j] - t.star[i]| are generated by the gp.predict
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
    gp.out <- predict.gp(theta, acv.model, dat, t.star)
  } else {
    tau.kk <- matrix.tau(t.star, t.star)
    gp.out <- list(y = rep(theta[1], m), 
                   cov = acv.model(theta[-c(1,2)], tau.kk))
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

# -----------------------------------------------------------
# Make a 'snake' plot of a Gaussian Process model

plot.snake <- function(dat, 
                       sigma = 1, 
                       add = FALSE, 
                       col.snake = NULL,
                       col.line = NULL,
                       ...) {

  # -----------------------------------------------------------
  # Inputs: 
  #   dat    - 3 column data frame containing
  #     t    - vector of n observation times
  #     y    - vector of n observations
  #     add  - (TRUE/FALSE) if TRUE then make a new plot
  #     dy   - vector of n 'errors' on observations
  # col.snake - (optional) colour for the 'snake'
  # col.line -(optional) colour for the line
  #
  # Value:
  #     none
  #
  # Description:
  #  Produces a 'snake' plot. This is a plot showing the mean of a process
  #  against time (as a thick line) and the variance of the process illustrated
  #  as a band. The lines goes through y(t) and the band covers y(t) +/-
  #  sigma*dy(t) where dy(t) is the standard deviation - this is actually
  #  sqrt(diag(covariance_matrix)). Sigma is a constant (>0) to allow one to
  #  plot e.g. the +/-2 std.dev band.
  #
  # History:
  #  06/07/16 - v0.1 - First working version
  #
  # Simon Vaughan, University of Leicester
  # -----------------------------------------------------------
  
  # check arguments
  if (!exists('dat')) {stop('** Missing dat.')}
  if (!("t" %in% names(dat))) {stop('** Missing dat$t.')}
  if (!("y" %in% names(dat))) {stop('** Missing dat$y.')}
  if (!("dy" %in% names(dat))) {
    if ("cov" %in% names(dat)) {
      dat$dy <- sqrt( abs( diag(dat$cov) ) )
    } else {
      stop('** Missing dat$dy and dat$cov - must include one of these.')
    }
  }
  
  # extract the data {t, y, dy} to plot
  t <- dat$t
  y <- dat$y
  dy <- dat$dy

  # define colours for the plot
  if (is.null(col.snake)) {
    col.snake <- rgb(255, 192, 203, 100, maxColorValue = 255)
  }
  if (is.null(col.line)) {
    col.line <- rgb(255, 50, 50, 255, maxColorValue = 255)
  }

  # prepare the snake (moving right along the bottom edge, then left along the top edge)
  xx <- c(t, rev(t))
  yy <- c(y - sigma*dy, rev(y + sigma*dy))
  
  # open a new plot if needed
  if (add != TRUE) {
    plot(xx, yy, type = "n", bty = "n", ...)
  }

  # now add the snake
  polygon(xx, yy, col = col.snake, border = NA)
  
  # now add the line through the centre
  lines(gp$t, gp$y, col = col.line, lwd=3)
  
  # all done
}

# -----------------------------------------------------------
