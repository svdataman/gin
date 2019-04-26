# -----------------------------------------------------------
#
# History:
#  25/04/19 - v0.1 - First working version
#
# Simon Vaughan, University of Leicester
# -----------------------------------------------------------

#' Load time series data into a standard format array.
#'
#' \code{ts_load} returns a standard-format data array
#'
#' @param filename  (string) name of data file
#' @param header    (logical) if true then line 1 is a header line.
#' @param csv       (logical) does the file use commas to seperate columns?
#' @param cols      (vector) list of columns to read
#'
#' @return
#'  array with three named columns \code{t, y, dy}.
#'
#' @section Notes:
#'  Load time series data from a multi-column ASCII file. The input data file
#'  should contain columns for \code{t}, \code{y} and \code{dy}.
#'  Here, \code{t} is the time variable, \code{y}
#'  is the variable measured and \code{dy} is its uncertainty (1-sigma error
#'  bar). If \code{dy} is not provided then it is assumed to be 0 for all data
#'  points. The output is a data frame of length N with columnes for \code{t},
#'  \code{y} and \code{dy}.
#'  If \code{csv} is set to \code{TRUE} then the file will be read as a CSV
#'  file. If \code{header} is set to \code{TRUE} then the first line of the
#'  file is assumed to be a header (and is then ignored).
#'  The optional input \code{cols} can specify which columns contain the
#'  \code{t}, \code{y} and \code{dy} data, e.g. \code{cols = c(2,4,5)}.
#'
#' @export
ts_load <- function(filename = NA, header = FALSE, csv = FALSE, cols = NA) {

  # check the data file exists
  if (is.na(filename)) return(NULL)
  if (!file.exists(filename)) {
    cat('** ts_load: File not found', fill = TRUE)
    return(NULL)
  }

  # load the data file
  if (csv == FALSE) {
    dat <- read.table(filename, header = header)
  } else {
    dat <- read.csv(filename, header = header)
  }

  # keep only specified columns
  if (!is.na(cols[1])) {
    dat <- dat[, cols]
  }

  # length of dataset
  N <- NROW(dat)
  if (N < 1) warning("File containts no data.")

  # number of columns
  M <- NCOL(dat)
  if (M > 3) dat <- dat[, 1:3]

  M <- NCOL(dat)
  if (M < 2) warning("File has <2 columns, not enough data.")

  # include missing error bars (as zeroes) if not provided.
  if (M == 2) {
    dy <- rep(0, N)
    dat <- cbind(dat, dy)
  }

  # name columns
  cnames <- c("t", "y", "dy")
  colnames(dat) <- cnames

  # return the data array in matrix format
  return(data.frame(dat))

}
# -----------------------------------------------------------
