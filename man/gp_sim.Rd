% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gp_functions.R
\name{gp_sim}
\alias{gp_sim}
\title{Generate a random GP realisation.}
\usage{
gp_sim(theta, acv.model = NULL, dat = NULL, t.star = NULL,
  N.sim = 1, plot = FALSE)
}
\arguments{
\item{theta}{(vector) parameters for covariance function
the first element is the mean value mu}

\item{acv.model}{(name) name of the function to compute ACV(tau | theta)}

\item{dat}{(matrix) an N * 3 matrix of data: 3 columns}

\item{t.star}{(vector) times at which to compute the simulation(s).}

\item{N.sim}{(integer) number of random vectors to draw.}

\item{plot}{(logical) add lines to an existing plot showing the results.}
}
\value{
An \code{n * N.sim} matrix, each column is a simulation, a pseudo-random
 \code{n}-length vector drawn from the Gaussian Process.
}
\description{
\code{gp_sim} returns one or more data vectors drawn from a Gaussian Process.
}
\section{Notes}{

  Simulate random realisations of a Gaussian Process (GP) with
  covariance defined by parameters \code{theta}, at times \code{t.star}.
  This is essentially eqn 2.22 of Rasmussen & Williams (2006):

  y.sim ~ N(mean = y.star, cov = C)

  Uses the function \code{rmvnorm} from the \code{mvtnorm} package to
  generate \code{n}-dimensional Gaussian vectors.

  If \code{dat} is supplied (a data frame/list containing \code{t} and
  \code{y}) then the conditional mean values at times \code{t.star},
  \code{y[t.star]} and covariance matrix \code{C} for lags \code{tau[i,j] =
  |t.star[j] - t.star[i]|} are generated by the \code{gp_conditional}
  function. See \code{gp_logLikelihood} for more about the input data form.

  If \code{dat} is not supplied we sample from the 'prior' GP, i.e. assume
  mean = \code{theta[1]} and covariance matrix given by \code{acv.model} and
  \code{theta} parameters.

  This function essentially just computes a covariance matrix given an
  autocovariance function (and its parameters), which may be conditional
  on some data, and uses this to draw a random Gaussian vector with
  \code{rmvnorm}.
}

\seealso{
\code{\link{gp_conditional}}, \code{\link{rmvnorm}}
}
