#' Main function of mapprox() in R
#'
#' \code{mapprox} is an wrapper of \code{linear_interpolation} and
#' \code{linear_interpolation_cpp}.
#'
#' This is the main function of \code{mapprox(..., use_cpp = FALSE)}.
#' Provided data are cleaned beforehand in \code{mapprox} and
#' then passed to this function.
#' \code{x} must be completely gridded data,
#' and duplicated values are not allowed.
#' Also, no NA in \code{x}, \code{y}, and \code{xout} is allowed.
#'
#' @param x Reference data of explanatory variables.
#'   This must be a data frame in which elements are all numeric.
#' @param y Reference data of a target variable.
#'   This must be a numeric vector with length \code{nrow(x)}.
#' @param xout Values of explanatory variables where interpolation is to
#'   take place. This must be a data frame and
#'   \code{ncol(xout)} must be the same as \code{ncol(x)}.
#' @param rule A scalar integer determining how interpolation outside the
#'   given explanatory variables works. This must be eigher
#'   \code{1}, \code{2}, or \code{3}.
#'   If \code{rule = 1}, \code{NA}s are returned for such points.
#'   If \code{rule = 2}, the values at the closest points in the region
#'   given as the reference data are returned.
#'   If \code{rule = 3}, linear extension of the interpolated result is applied.
#' @param verbose Logical. Should progress of the process printed?
#'
#' @return A numeric vector of interpolated result.
#'
#' @seealso [mapprox()]

linear_interpolation <- function(x, y, xout, rule, verbose) {

  # Interpolation (1): Preparation ---------------------------------------------
  t1 <- proc.time()[3]
  if (verbose) cat('Step 3/4: Preparation...     ')

  # Converting x to a list and y to an array
  x_list <- lapply(x, function(x0) sort(unique(x0)))
  y_arr <- tapply(X = y, INDEX = as.list(x), FUN = identity)
  names(dimnames(y_arr)) <- colnames(x)

  # When rule = 2, values of xout is replaced with closest value of x
  # if they are out of range of x.
  if (rule == 2) {
    xout <- as.data.frame(mapply(
      function(x0, xout0) pmax(pmin(xout0, max(x0)), min(x0)),
      x0 = x_list, xout0 = xout))
  }

  # Where does each data of xout place? Which grid expanded by x?
  # An element of `smaller_indices` (i_rc in row r and column c) indicate
  # that the xout[r, c] is between i_rc_th and (i_rc + 1)_th values of x
  # (x_list[[c]][i_rc + (0:1)]).
  smaller_indices <- mapply(function(x, vals) {

    i <- sapply(vals, function(val) sum(val - x >= 0)) # faster
    # i <- sapply(vals, function(val) which(val - x < 0)[1] - 1) # slower

    # When rule = 3, xout placed out of the region of x are treated
    # as if they are in the closest grid.
    pmax(pmin(i, length(x) - 1), 1)

  }, x = x_list, vals = xout)

  # Lengths of grids where data points of xout placed
  region_length <- mapply(
    function(x, i) x[i + 1] - x[i],
    x = x_list, i = as.data.frame(smaller_indices))

  t2 <- proc.time()[3]
  if (verbose) cat(' Done. (', round(t2 - t1, 1), 's)\n', sep = '')
  t1 <- t2


  # Interpolation (2): calculation of yout -------------------------------------
  if (verbose) cat('Step 4/4: Calculating yout...')

  # All 0/1 patterns corresponding to each explanatory variable
  added_indices <- do.call(expand.grid, rep(list(0:1), ncol(x)))
  colnames(added_indices) <- colnames(x)

  # Loop for each grid points (x) surrounding xout
  yout <- rep(0, nrow(xout))
  print_progress <- FALSE
  n_loop <- nrow(added_indices)
  for (j in seq(n_loop)) {
    # j <- 1

    added_index_j <- unlist(added_indices[j, ])

    # y of grid point j
    y_indices_j <- smaller_indices + added_index_j[col(smaller_indices)]
    y_j <- y_arr[y_indices_j]

    # components of weights for grid point j
    rev_added_index_j <- 1 - added_index_j
    x_rev_indices_j <- smaller_indices + rev_added_index_j[col(smaller_indices)]
    sign_adjuster_j <- added_index_j * 2 - 1
    weight_raw_j <- mapply(
      FUN = function(x, reglen, ref_x, x_rev_ind, sign_adj) {
        sign_adj * (x - ref_x[x_rev_ind]) / reglen
      },
      x         = xout,
      reglen    = as.data.frame(region_length),
      ref_x     = x_list,
      x_rev_ind = as.data.frame(x_rev_indices_j),
      sign_adj  = sign_adjuster_j)

    # weights for grid point j (slower)
    # weight_j <- apply(X = weight_raw_j, MARGIN = 1, FUN = prod)

    # weights for grid point j (faster)
    weight_j <- rep(0, nrow(weight_raw_j))
    is_0 <- rowSums(weight_raw_j == 0) > 0
    weight_raw_not0_j <- weight_raw_j[!is_0, , drop = FALSE]
    prodabs_weight_raw_j <- exp(rowSums(log(abs(weight_raw_not0_j))))
    n_minus_weight_raw_j <- rowSums(sign(weight_raw_not0_j) == -1)
    weight_j[!is_0] <- (-1) ^ (n_minus_weight_raw_j %% 2) * prodabs_weight_raw_j

    # Renew the answer
    yout <- yout + y_j * weight_j

    if (!print_progress && verbose && (proc.time()[3] - t1 > 10)) {
      print_progress <- TRUE
      cat('\n')
    }

    if (print_progress) {
      finish_time <- Sys.time() + (proc.time()[3] - t1) / j * (n_loop - j)
      cat('\r  Loop ', j, '/', n_loop, ', finish time prediction: ',
          as.character(round(finish_time)), sep = '')
    }

  }
  if (print_progress) cat('\n                             ')

  # If rule = 1, the returned values should be NA
  # for xout placed out of range of x.
  if (rule == 1) {
    is_out_of_range <- mapply(
      function(x0, xout0) xout0 < min(x0) | xout0 > max(x0),
      x0 = x_list, xout0 = xout)
    yout[rowSums(is_out_of_range) > 0] <- NA
  }

  t2 <- proc.time()[3]
  if (verbose) cat(' Done. (', round(t2 - t1, 1), 's)\n', sep = '')
  names(yout) <- NULL
  dim(yout) <- NULL
  yout

}

#' Linear interpolation of a multivariable function.
#'
#' \code{mapprox} is an extension of \code{approx} to deal with multivariable x.
#'
#' @param x Reference data of explanatory variables.
#'   This must be either a matrix or data frame
#'   in which elements are all numeric.
#'   A vector is also allowed in an univariable case.
#' @param y Reference data of a target variable.
#'   This must be a numeric vector with length \code{nrow(x)}.
#' @param xout Values of explanatory variables where interpolation is to
#'   take place. This must be a matrix or data frame and
#'   \code{ncol(xout)} must be the same as \code{ncol(x)}.
#'   As for \code{x}, a vector is also allowed in an univariable case.
#' @param rule A scalar integer determining how interpolation outside the
#'   given explanatory variables works. This must be eigher
#'   \code{1}, \code{2}, or \code{3}.
#'   If \code{rule = 1}, \code{NA}s are returned for such points.
#'   If \code{rule = 2}, the values at the closest points in the region
#'   given as the reference data are returned.
#'   If \code{rule = 3}, linear extension of the interpolated result is applied.
#' @param use_cpp Logical. Do you use C++ for faster calculation?
#' @param verbose Logical. Should progress of the process printed?
#'
#' @return A list with components \code{x}, \code{y}, \code{yout}.
#'   The elements \code{x} and \code{y} are identical to
#'   input \code{x} and \code{y}, respectively,
#'   if all data are allowed as appropriate input.
#'   Otherwise, they are subsets of input reference data
#'   without irregular values (NA, non-gridded values, and duplicated values).
#'   The element \code{yout} is a numeric vector with length \code{nrow(xout)}
#'   containing interpolated results.
#'
#' @importFrom Rcpp sourceCpp
#' @useDynLib mapprox, .registration = TRUE
#' @seealso [approx()]
#'
#' #' @examples
#' # Example 1. A one-variable case: comparison with approx().
#' x1 <- 1:10
#' y1 <- c(1:3, 8:10, 7:4)
#' xout1 <- seq(0, 11, 0.2)
#' res1_approx <- approx(x1, y1, xout1, method = 'linear', rule = 2)
#' res1_mapprox <- mapprox(x1, y1, xout1, rule = 2)
#' plot(res1_approx, pch = 16, cex = 0.5, xlab = 'x', ylab = 'y')
#' points(xout1, res1_mapprox$yout)
#' points(x1, y1, pch = 16, cex = 2, col = rgb(1, 0, 0, 0.2))
#' legend('topleft', legend = c('Data', 'approx', 'mapprox'),
#'        col = c(rgb(1, 0, 0, 0.2), 'black', 'black'),
#'        pch = c(16, 16, 1), pt.cex = c(2, 1, 1))
#'
#' # Example 2. A two-variable case.
#' x2 <- expand.grid(V1 = 11:15, V2 = 101:103)
#' y2 <- with(x2, (V1 - 11.5) ^ 2 + V2 * 5)
#' plot(x2$V1, y2, col = x2$V2 - 100, pch = 16, xlab = 'V1',
#'      xlim = c(11, 16), ylim = c(min(y2), 530))
#' legend('topleft', legend = paste0('V2=', 101:103),
#'        col = 1:3, pch = 16, ncol = 3)
#' xout2 <- expand.grid(V1 = seq(11, 16, 0.1), V2 = c(101.5, 102.8))
#' res2 <- mapprox(x2, y2, xout2, rule = 2)
#' points(xout2$V1, res2$yout, col = ifelse(xout2$V2 == 101.5, 4, 5),
#'        pch = 16, cex = 0.4)
#' legend('bottomright', legend = paste0('V2=', c(101.5, 102.8)),
#'        col = 4:5, pch = 16, ncol = 2)
#'
#' # Example 3. Treatment of NA, non-gridded values, and duplicated values.
#' x3 <- expand.grid(V1 = 11:15, V2 = 101:103) # correct data format
#' y3 <- with(x3, (V1 - 12) ^ 2 + V2 * 5)
#' x3 <- rbind(x3,
#'             data.frame(V1 = NA, V2 = 103), # NA in x
#'             data.frame(V1 = 16, V2 = 103), # non-gridded value
#'             data.frame(V1 = 14, V2 = 103), # duplicated value
#'             data.frame(V1 = 12, V2 = 103)) # NA in y
#' y3 <- c(y3, 500, 531, 525, NA)
#' plot(x3$V1, y3, col = x3$V2 - 100, pch = 16, xlab = 'V1')
#' legend('topleft', legend = paste0('V2=', 101:103),
#'        col = 1:3, pch = 16, ncol = 3)
#' xout3 <- expand.grid(V1 = seq(11, 16, 0.1), V2 = 101:103)
#' res3 <- mapprox(x3, y3, xout3, rule = 2) # Three warnings
#' res3$x
#' res3$y
#' points(xout3$V1, res3$yout, col = xout3$V2 - 100, pch = 16, cex = 0.4)
#' # NA and non-gridded data are ignored and duplicated data are averaged
#'
#' # Example 4. Comparison of rules.
#' x4 <- expand.grid(V1 = 11:14, V2 = 101:102)
#' y4 <- with(x4, (V1 - 11.5) ^ 2 + V2 * 5)
#' xout4 <- expand.grid(V1 = seq(11, 15, 0.1), V2 = 101.5)
#' plot(x4$V1, y4, col = x4$V2 - 100, pch = 16, xlab = 'V1',
#'      xlim = range(xout4$V1))
#' legend('topleft', legend = paste0('V2=', 101:102),
#'        col = 1:2, pch = 16, ncol = 2)
#' res4_1 <- mapprox(x4, y4, xout4, rule = 1)
#' res4_2 <- mapprox(x4, y4, xout4, rule = 2)
#' res4_3 <- mapprox(x4, y4, xout4, rule = 3)
#' points(xout4$V1, res4_1$yout, col = 3, pch = 3, cex = 2)
#' points(xout4$V1, res4_2$yout, col = 4, pch = 15)
#' points(xout4$V1, res4_3$yout, col = 5, pch = 16)
#' legend('bottomright', legend = paste0('rule=', 1:3),
#'        col = 3:5, pch = c(3, 15, 16), ncol = 3)

mapprox <- function(x, y, xout, rule = 1, use_cpp = FALSE, verbose = FALSE) {

  # Checking arguments ---------------------------------------------------------
  t1 <- proc.time()[3]
  n_step <- if (use_cpp) 3 else 4
  if (verbose) cat('Step 1/', n_step, ': Checking inputs... ', sep = '')

  # Checking `x`
  if (!(is.vector(x) || is.matrix(x) || is.data.frame(x))) {
    stop('x must be a vector, matrix, or data frame.')
  }
  if (is.vector(x)) x <- as.matrix(x, ncol = 1)
  if (is.matrix(x)) x <- as.data.frame(x)
  if (any(!sapply(x, is.numeric))) stop('x contains values not numeric.')

  # Checking `y`
  if (!is.vector(y)) stop('y must be a vector.')
  if (!is.numeric(y)) stop('y must be numeric.')
  if (length(y) != nrow(x)) {
    stop('length(y) must be the same as the number of variables in x.')
  }

  # checking `xout`
  if (!(is.vector(xout) || is.matrix(xout) || is.data.frame(xout))) {
    stop('xout must be a vector, matrixout, or data frame.')
  }
  if (is.vector(xout)) xout <- as.matrix(xout, ncol = 1)
  if (is.matrix(xout)) xout <- as.data.frame(xout)
  if (any(!sapply(xout, is.numeric))) stop('xout contains values not numeric.')
  if (ncol(xout) != ncol(x)) {
    stop('The number of explanatory variables must be the same for x and xout.')
  }

  # checking `rule`
  if (!(identical(rule, 1) || identical(rule, 2) || identical(rule, 3))) {
    stop('rule must be either 1, 2, or 3.')
  }

  # checking `use_cpp`
  if (!is.logical(use_cpp)) stop('use_cpp must be logical.')
  if (length(use_cpp) != 1) stop('use_cpp must be length 1.')

  # checking `verbose`
  if (!is.logical(verbose)) stop('verbose must be logical.')
  if (length(verbose) != 1) stop('verbose must be length 1.')

  t2 <- proc.time()[3]
  if (verbose) cat(' Done. (', round(t2 - t1, 1), 's)\n', sep = '')
  t1 <- t2


  # Reshaping data -------------------------------------------------------------
  if (verbose) cat('Step 2/', n_step, ': Reshaping data...  ', sep = '')

  # Omitting NA
  is_na_x <- rowSums(is.na(x)) > 0
  is_na_y <- is.na(y)
  is_na <- is_na_x | is_na_y
  x <- x[!is_na, , drop = FALSE]
  y <- y[!is_na]
  if (any(is_na)) {
    warning(
      'There were ', sum(is_na_x), ' NAs in x and ',
      sum(is_na_y), ' NAs in y. ',
      sum(is_na), ' data were omitted in total.')
  }

  # Checking if x is a gridded data
  if (ncol(x) > 1) {
    is_gridded <- extract_gridded_x(x)
    x <- x[is_gridded, , drop = FALSE]
    y <- y[is_gridded]
    if (any(!is_gridded)) {
      warning(sum(!is_gridded), ' values of x were not gridded and so dropped.')
    }
  }

  # Averaging y with duplicated x
  duplication_id <- detect_duplicated_x(x)
  if (any(!is.na(duplication_id))) {
    warning(
      sum(!is.na(duplication_id)),
      ' values of x were duplicated and their y values are averaged.'
    )
  }
  xy <- remove_duplicated_data(x, y, duplication_id)
  x <- xy$x
  y <- xy$y

  # Stopping the process if no data remained after cleaning
  if (nrow(x) == 0) stop('No data remained after data cleaning.')

  t2 <- proc.time()[3]
  if (verbose) cat(' Done. (', round(t2 - t1, 1), 's)\n', sep = '')
  t1 <- t2


  # Interpolation --------------------------------------------------------------
  # Removing xout with NA
  is_na_xout <- rowSums(is.na(xout)) > 0
  if (all(is_na_xout)) stop('All points of xout contain NA.')
  xout_nona <- xout[!is_na_xout, , drop = FALSE]
  yout <- rep(NA, nrow(xout))

  # Interpolation
  ord <- do.call(order, as.list(x))
  x <- x[ord, , drop = FALSE]
  y <- y[ord]
  if (use_cpp) {

    if (verbose) cat('Step 3/', n_step, ': Interpolation ...  ', sep = '')
    yout[!is_na_xout] <- linear_interpolation_cpp(x, y, xout_nona, rule, verbose)
    t2 <- proc.time()[3]
    if (verbose) cat(' Done. (', round(t2 - t1, 1), 's)\n', sep = '')

  } else {

    yout[!is_na_xout] <- linear_interpolation(x, y, xout_nona, rule, verbose)

  }

  list(x = x, y = y, yout = yout)

}

