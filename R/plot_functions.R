# -----------------------------------------------------------
# History:
#  06/07/16 - v0.1 - First working version
#  12/07/16 - v0.2 - Improved use of colour inputs
#
# Simon Vaughan, University of Leicester
# -----------------------------------------------------------

#' Make a 'snake' plot of a Gaussian Process model
#'
#' \code{plot_snake} plots the mean and N-sigma band for a GP.
#'
#' Produces a 'snake' plot. This is a plot showing the mean of a Gaussian
#'  process against time (as a thick line) and the variance of the process
#'  illustrated as a band. The lines goes through mu(t) and the band covers
#'  mu(t) +/- sigma*dy(t) where dy(t) is the standard deviation - this is
#'  actually sqrt(diag(covariance_matrix)). Sigma is a constant (>0) to allow
#'  one to plot e.g. the +/-2 std.dev band.
#'
#' @param dat (data frame) data for time series.
#' @param sigma (float) Plot the +/-sigma band (>0).
#' @param add (logical) if TRUE then start a new plot
#' @param col.line (colour) Colour of the mean line.
#' @param col.fill (colour) Colour of the +/-sigma band.
#' @param col.border (colour) Colour of +/-sigma limits.
#' @param ... (anything) any other graphical keywords to be passed to
#'                        \code{plot(...)}
#'
#' @return
#'   None.
#'
#' @seealso \code{\link{plot_ts}}
#'
#' @export
plot_snake <- function(dat,
                       sigma = c(1),
                       add = FALSE,
                       col.line = 1,
                       col.fill = NULL,
                       col.border = NULL,
                       ...) {


  # check arguments
  if (!exists('dat')) {stop('** Missing dat.')}

  # is this a gp object (list format) or a data array?
  is.gp <- mode(dat) == "list"

  # if data array, check the contents
  if (is.gp == FALSE) {
    if (!("t" %in% colnames(dat))) {stop('** Missing dat$t.')}
    if (!("y" %in% colnames(dat))) {stop('** Missing dat$y.')}
    if (!("dy" %in% colnames(dat))) {stop('** Missing dat$dy.')}
  } else {
    if (!("t" %in% names(dat))) {stop('** Missing dat$t.')}
    if (!("y" %in% names(dat))) {stop('** Missing dat$y.')}
    if (!("dy" %in% names(dat))) {
      if ("cov" %in% names(dat)) {
        dat$dy <- sqrt( abs( diag(dat$cov) ) )
      } else {
        stop('** Missing dat$dy and dat$cov - must include one of these.')
      }
    }
  }

  # extract the data {t, y, dy} to plot
  t <- dat$t
  y <- dat$y
  dy <- dat$dy

  # define colours for the plot
  cols <- col.line
  if (mode(col.line) == "numeric") {
    col.line <- RColorBrewer::brewer.pal(8, "Set1")[cols]
  }

  if (is.null(col.fill)) {
    col.fill <- col.line
  }
  if (is.null(col.border)) {
    col.border <- NA
  }
  col.fill <- rgb(t(col2rgb(col.fill)), alpha = 50,
                  maxColorValue = 255)

  # loop over multiple sigma bands, if needed
  n.bands <- length(sigma)
  for (i in 1:n.bands) {
    # prepare the snake (moving right along the bottom edge, then left along the
    # top edge)
    xx <- c(t, rev(t))
    yy <- c(y - sigma[i] * dy, rev(y + sigma[i] * dy))

    # open a new plot if needed
    if (add != TRUE) {
      if (i == 1)
        plot(xx, yy, type = "n", bty = "n", ...)
    }

    # now add the snake
    polygon(xx,
            yy,
            col = col.fill,
            border = col.border,
            lwd = 2)
  }

  # now add the line through the centre
  lines(t, y, col = col.line, lwd = 3)

  # all done
}

# -----------------------------------------------------------

# -----------------------------------------------------------
#' Nice plot of time series data
#'
#' \code{plot_ts} plots uneven time series data frame.
#'
#' Makes a nice plot of some time series data. The data are input as a data
#' frame
#'
#' @param ts1 (data frame) data for time series.
#' @param ts2 (data frame) option data for second time series.
#' @param xlim,ylim (vectors) numeric vectors of length 2, giving the x and y
#'                     coordinates ranges.
#' @param xlab,ylab (strings) titles for the x, y axes.
#' @param grid (logical) plot a grid?
#' @param ... (anything) any other graphical keywords to be passed to
#'                        \code{plot(...)}
#' @param cex (float) A numerical value giving the amount by which plotting
#'                    text and symbols should be magnified relative to the
#'                    default.
#' @param col (colour) Specification of the point/line colour.
#' @param pch (integer) An integer specifying the symbol
#'                      to be used as the default in plotting points.
#' @param main (string) An overall title
#' @param type (character) What type of plot. See \code{?plot}.
#'             Options are "p", "b", "l", "n".
#' @param lwd (float) line width.
#' @param extend (float) Extend the time axis by this factor.
#' @param grid.lwd (float) The line width (>0) for the grid.
#' @param grid.col (colour) The colour of the grid.
#' @param horizon (logical or float) make a 'horizon' style plot?
#' @param error (logical) plot errors if available?
#' @param cex.lab (float) The magnification to be used for x and y labels
#'                 relative to the current setting of cex.
#' @inheritParams plot_snake
#'
#' @return
#'   None.
#'
#' @section Notes:
#' Input data frames: note that the input data \code{dat} is not a traditional
#' \code{R} time series object. Such objects are only suitable for regularly
#' sampled data. These plotting functions are designed to work with data of
#' arbitrary sampling, we therefore need to explicitly list times and values.
#' The input objects are therefore data.frames with at least two columns which
#' much be called \code{t} (time) and \code{y} (value). An error of the value
#' may be provided by a \code{dy} column. If the times are not points but bins
#' of finite width, the width of each time bin may be specified with the
#' \code{dt} column. Any other columns are ignored.
#'
#' @seealso \code{\link{plot_ts}}
#'
#' @examples
#' plot_ts(drw)
#'
#' @export
plot_ts <- function(ts1, ts2,
                    xlim = NULL,
                    ylim = NULL,
                    xlab = "time",
                    ylab = "y(t)",
                    grid = TRUE,
                    add = FALSE,
                    pch = 16,
                    cex = 1,
                    lwd = 3,
                    col = NULL,
                    main = "",
                    extend = 1.0,
                    error = TRUE,
                    type = "p",
                    grid.lwd = 1,
                    grid.col = "lightgray",
                    horizon = FALSE,
                    split = FALSE,
                    cex.lab = 1.5,
                     ...) {

  # check arguments
  if (!exists('ts1')) {stop('** Missing ts1.')}
  if (!("y" %in% colnames(ts1))) stop('** Missing ts1$y.')
  n <- length(ts1$y)
  y <- ts1$y
  t <- 1:n

  if ("t" %in% names(ts1)) t <- ts1$t

  dy <- array(0, dim = n)
  dt <- array(0, dim = n)
  if ("dy" %in% colnames(ts1)) dy <- ts1$dy
  if ("dt" %in% colnames(ts1)) dt <- ts1$dt
  if (error == FALSE) {
    dy <- 0
    dt <- 0
  }

  # set ranges
  if (is.null(xlim)) {
    xlim <- range(t)
    if (!missing(ts2)) xlim <- range(c(xlim, range(ts2$t)))
  }
  trange <- diff(xlim)
  xlim <- xlim + c(-1, 1) * (extend-1)*trange/2

  if (mode(col) == "character") {
    cols <- col
  } else {
    cols <- RColorBrewer::brewer.pal(8, "Set1")
  }
  icol <- 1
  if (is.numeric(col)) icol <- col

  # store the current graphical parameters
  retire <- par(no.readonly = TRUE)

  m <- 1
  if (!missing(ts2)) m <- 2

  if (add != TRUE) {
    if (split == TRUE) {
      par(mfrow = c(m, 1))
#      par(mar = c(0, 0, 0, 0), oma = c(5, 5, 0.5, 0.5))
    }
  }

  yfix <- TRUE
  if (is.null(ylim)) yfix <- FALSE
  if (yfix != TRUE) {
    if (!missing(ts2) & split != TRUE) {
      ylim <- range(c(ts1$y, ts2$y))
      yfix <- TRUE
    } else {
      ylim <- range(ts1$y)
    }
  }

  for (i in 1:m) {

    if (i > 1) {
      y <- ts2$y
      n <- length(y)
      t <- 1:n
      if ("t" %in% colnames(ts2)) t <- ts2$t
      dy <- array(0, dim = n)
      if ("dy" %in% colnames(ts2)) dy <- ts2$dy
      dt <- array(0, dim = n)
      if ("dt" %in% colnames(ts2)) dt <- ts2$dt
      icol <- icol + 1
    }

    # open a new plot if needed
    if (add != TRUE) {
      if (i == 1 | split == TRUE) {
        if (yfix != TRUE)
          ylim <- range(y)
        plot(0, 0, xlim = xlim, ylim = ylim, axes = FALSE,
             type = "n", bty = "n", xlab = "", ylab = ylab,
             main = main, cex.lab = cex.lab, ...)
        axis(2, lty = 0, las = 1)
        abline(h = par("usr")[3])
        if (grid == TRUE) grid(lwd = grid.lwd, col = grid.col)
      }
    }


  # draw 'horizon' if requested
  bottom <- NULL
  if (mode(horizon) == "logical" & horizon == TRUE)
    bottom <- par("usr")[3]
  if (mode(horizon) == "numeric")
    bottom <- horizon

  col.fill <- rgb(t(col2rgb(cols[icol])), alpha = 50,
                  maxColorValue = 255)

  if (!is.null(bottom)) {
    pt <- c(t[1], t, t[n])
    py <- c(bottom, y, bottom)
    polygon(pt, py, col = col.fill, border = NA)
  }

  # draw points or lines
    if (type == "p" | type == "b") {
      points(t, y, pch = pch, cex = cex, col = cols[icol])
    }
    if (type == "l" | type == "b") {
      lines(t, y, col = cols[icol], lwd = lwd)
    }

  # include errors [optionally]
    if (error == TRUE) {
      mask <- which(dy > 0)
      segments(t[mask], y[mask]-dy[mask], t[mask], y[mask]+dy[mask],
               col = cols[icol], lwd = max(1, lwd-1))
      mask <- which(dt > 0)
      segments(t[mask]-dt[mask]/2, y[mask], t[mask]+dt[mask]/2, y[mask],
               col = cols[icol], lwd = max(1, lwd-1))
    }
  }
  axis(1)
  mtext(xlab, side = 1, line = 3, cex = cex.lab)

  # restore graphical parameters
  # par(retire)

}
