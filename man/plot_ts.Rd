% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_functions.R
\name{plot_ts}
\alias{plot_ts}
\title{Nice plot of time series data}
\usage{
plot_ts(ts1, ts2 = NULL, xlim = NULL, ylim = NULL, xlab = "time",
  ylab = "y(t)", grid = TRUE, add = FALSE, pch = 16, cex = 1,
  lwd = 3, col = NULL, main = "", extend = 1, error = TRUE,
  type = "p", grid.lwd = 1, grid.col = "lightgray",
  horizon = FALSE, split = FALSE, cex.lab = 1.5, ...)
}
\arguments{
\item{ts1}{(matrix) data for time series.}

\item{ts2}{(data frame) option data for second time series.}

\item{xlim, ylim}{(vectors) numeric vectors of length 2, giving the x and y
coordinates ranges.}

\item{xlab, ylab}{(strings) titles for the x, y axes.}

\item{grid}{(logical) plot a grid?}

\item{add}{(logical) if TRUE then start a new plot}

\item{pch}{(integer) An integer specifying the symbol
to be used as the default in plotting points.}

\item{cex}{(float) A numerical value giving the amount by which plotting
text and symbols should be magnified relative to the
default.}

\item{lwd}{(float) line width.}

\item{col}{(colour) Specification of the point/line colour.}

\item{main}{(string) An overall title}

\item{extend}{(float) Extend the time axis by this factor.}

\item{error}{(logical) plot errors if available?}

\item{type}{(character) What type of plot. See \code{?plot}.
Options are "p", "b", "l", "n".}

\item{grid.lwd}{(float) The line width (>0) for the grid.}

\item{grid.col}{(colour) The colour of the grid.}

\item{horizon}{(logical or float) make a 'horizon' style plot?}

\item{cex.lab}{(float) The magnification to be used for x and y labels
relative to the current setting of cex.}

\item{...}{(anything) any other graphical keywords to be passed to
\code{plot(...)}}
}
\value{
None.
}
\description{
\code{plot_ts} plots uneven time series data frame.
}
\details{
Makes a nice plot of some time series data. The data are input as a data
frame
}
\section{Notes}{

Input data: note that the input data \code{dat} is not a traditional
\code{R} time series object. Such objects are only suitable for regularly
sampled data. These plotting functions are designed to work with data of
arbitrary sampling, we therefore need to explicitly list times and values.
The input objects is a matrix with at least two columns
called \code{t} (time) and \code{y} (value). An error of the value
may be provided by a \code{dy} column. If the times are not points but bins
of finite width, the width of each time bin may be specified with a
\code{dt} column. Any other columns are ignored.
}

\examples{
plot_ts(drw)

}
\seealso{
\code{\link{plot_ts}}
}
