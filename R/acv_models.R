# -----------------------------------------------------------
# define power spectrum model in terms of parameters 'theta'

psd.model <- function(theta, f) {
  
  # -----------------------------------------------------------
  # psd.model
  # Inputs: 
  #   theta - vector of (hyper-)parameters for PSD
  #   f     - frequencies at which to compute PSD
  #
  # Value:
  #  psd   - power density at frequencies f
  #
  # Description:
  #  Compute the Power Spectral Density (PSD) at input frequencies
  # based on parameters theta. 
  #
  # History:
  #  16/12/13 - First working version
  #
  # Simon Vaughan, University of Leicester
  # -----------------------------------------------------------
  
  # check arguments
  if (missing(theta)) {stop('** Missing THETA argument in PSD.MODEL')}
  if (missing(f)) {stop('** Missing F argument in PSD.MODEL')}
  
  # PSD is symmetric   
  f <- abs(f)
  df <- f[2]-f[1]
  
  # parameters 
  a <- theta[2]
  n <- 10^theta[1] 
  
  f.min <- 1e-2
  
  # define power law model  
  psd <- n * f^(-a)
  mask <- (f < f.min)
  min.j <- which.min(f[f > f.min])
  psd[mask] <- (psd[f > f.min])[min.j]
  
  # ensure DC power = 0  
  #  ii <- which.min(f)
  #  if (f[ii] == 0) { psd[ii] <- 0 }
  
  # finished
  return(psd)
}

# -----------------------------------------------------------
# compute ACF from power spectral model

acf.model <- function(theta, lags, diagnose=FALSE) {
  
  # -----------------------------------------------------------
  # acf.model
  # Inputs: 
  #   theta - vector of (hyper-)parameters for ACF/PSD
  #   lags  - lags at which to compute ACF
  #
  # Value:
  #  acov   - Auto-covariance ACF(lags)
  #
  # Description:
  #  Compute the auto-covariance function (ACF) at input lags
  # based on parameters theta. In practice the ACF is the Fourier
  # transform of a PSD model. First we define the Fourier frequencies,
  # then compute the PSD model at these frequencies, then we include
  # the negatives frequencies - PSD(-f) = PSD(f) - and perform an FFT. 
  # This gives the ACF values. We keep only the ACF values at zero and 
  # positive lags
  #
  # History:
  #  16/12/13 - First working version
  #
  # Simon Vaughan, University of Leicester
  # -----------------------------------------------------------
  
  # check arguments
  if (missing(theta)) {stop('** Missing THETA argument in ACF.MODEL')}
  if (missing(lags)) {stop('** Missing LAGS argument in ACF.MODEL')}
  
  #   # compute corresponding frequencies
  #   # Note: we need to use twice the usual Fourier frequencies here
  #   n.t <- length(lags)
  #   duration <- lags[n.t]
  #   
  #   df <- 1/duration/2
  #   nf <- n.t          
  #   f.max <- nf * df
  #   f <- seq(0, f.max, by=df)  
  #   
  #   # compute PSD model
  #   psd <- psd.model(theta, f)
  #   
  #   # make sure the input PSD is defined for equally spaced positive 
  #   # frequencies 'f' running 0,1,...,N.f.
  #   n.f <- length(psd)
  #   
  #   # include the negative frequencies (PSD is symmetric)
  #   # factor of 1/2 needed when we move from one-sided to
  #   # two-sided PSD
  #   mask <- (1:(n.f-2))+1
  #   psd.sym <- c(psd, rev(psd[mask])) / 2
  #   
  #   # compute the FFT
  #   # XXX check normalisation term XXX
  #   dt <- 1/f.max/2
  #   acov <- Re(fft(psd.sym)) / length(psd.sym) / dt
  #   
  #   # keep only the positive lags (since symmetric)
  #   acov <- acov[1:(n.f-1)]
  #   
  #   if (diagnose == TRUE) {
  #     cat('-- theta: ', theta, fill=TRUE)
  #     cat('-- lags:  ', lags, fill=TRUE)
  #     cat('-- acf:   ', acov, fill=TRUE)
  #   }
  
  #  plot(psd)
  #  plot(acov)
  
  acov <- (10^theta[1])*exp(-lags/theta[2])
  
  return(acov)
  
}

# -----------------------------------------------------------
# build the covariance matrix from the PSD model
#
#  tau = [N.i, N.j] matrix of lags tau[i,j] = |t.i - t.j|

acv <- function(theta, tau) {
  
  # -----------------------------------------------------------
  # matrix.acf
  # Inputs: 
  #   theta - vector of (hyper-)parameters for ACF/PSD
  #   tau   - N*M array of lags at which to compute ACF
  #
  # Value:
  #  acf.ij - array of auto-covariances ACF(tau)
  #
  # Description:
  #  Compute the auto-covariance function (ACF) at input lags tau,
  # based on parameters theta. Tau is an N*M array of lags. We first
  # compute the 1-dimensional ACF(lags) and then use linear interpolation
  # to compute the ACF at every lag of the array tau[i,j] = |t1[j] - t2[i]|
  #
  # Note: we cannot simply assume tau is a cyclic, or square matrix,
  # hence we need to compute the ACF at each tau value explicitly.
  #
  # History:
  #  16/12/13 - First working version
  #
  # Simon Vaughan, University of Leicester
  # -----------------------------------------------------------
  
  # use package Matrix for Matrix classes and functions
  require(Matrix)
  
  # check arguments
  if (missing(theta)) {stop('** Missing THETA argument in MATRIX.ACF')}
  if (missing(tau)) {stop('** Missing TAU argument in MATRIX.ACF')}
  if (length(dim(tau)) != 2) {stop('** TAU is not 2d array in MATRIX.ACF')}
  
  # ranges for computing ACF
  t.range <- range(tau[tau > 0])
  dt.min <- t.range[1]
  dt.max <- t.range[2]
  n.i <- NROW(tau)
  n.j <- NCOL(tau)
  
  # vector of time delays for computing 1-dimensional ACF(dt)
  # XXX check this XXX
  n.step <- max(n.i, n.j, dt.max/dt.min+1)
  n.step <- n.step - n.step%%2             #  make even length
  lags <- seq(0, dt.max, length.out=n.step)
  
  # compute 1-dimensional ACF
  acf.i <- acf.model(theta, lags)
  
  # now interpolate this to give ACF at tau_ij
  # and reshape into an N.i*N.j matrix
  acf.func <- approxfun(lags, acf.i, yright=0)
  acf.ij <- acf.func(tau)
  dim(acf.ij) <- c(n.i, n.j)
  acf.ij <- Matrix(acf.ij)
  
  # all done
  return(acf.ij)
  
}
