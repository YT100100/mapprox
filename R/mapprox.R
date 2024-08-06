#' Linear interpolation of a multivariable function.
#'
#' \code{mapprox} is an extension of \code{approx} to deal with multivariable x.
#'
#' @param x Reference data of explanatory variables.
#'   This must be a vector, matrix or data frame in which elements are all numeric.
#' @param y Reference data of a target variable.
#'   This must be a numeric vector with length \code{nrow(x)}.
#' @param xout Values of explanatory variables where interpolation is to
#'   take place. This must be a matrix or data frame and
#'   \code{ncol(xout)} must be the same as \code{ncol(x)}.
#' @param rule A scalar integer determining how interpolation outside the
#'   given explanatory variables works. This must be eigher
#'   \code{1}, \code{2}, or \code{3}.
#'   If \code{rule = 1}, \code{NA}s are returned for such points.
#'   If \code{rule = 2}, the values at the closest points in the region
#'   given as the reference data are returned.
#'   If \code{rule = 3}, linear extension of the interpolated result is applied.
#'
#' @return A list with components \code{x} and \code{y}.
#'   The element \code{x} is identical to input \code{xout},
#'   if its all rows are allowed as appropriate input.
#'   The element \code{y} is a numeric vector with length \code{nrow(xout)},
#'   containing interpolated results.
#'
#' @seealso [approx()]
#'
#' @examples
#' # An example of one variable: comparison with approx().
#' x <- 1:10
#' y <- c(1:3, 8:10, 7:4)
#' xout <- seq(0, 11, 0.2)
#' res_approx <- approx(x, y, xout, method = 'linear', rule = 2)
#' res_mapprox <- mapprox(x, y, xout, rule = 2)
#' plot(res_approx, pch = 16, cex = 0.5)
#' points(res_mapprox)
#' points(x, y, pch = 1, cex = 2)
#' # legend()
#'
#'

mapprox <- function(x, y, xout, rule = 1) {

  # Example 1: One variable
  # x <- 1:10
  # y <- c(1:3, 8:10, 7:4)
  # xout <- seq(0, 11, 0.2)
  # rule <- 2


  # Checking arguments ---------------------------------------------------------
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


  # Reshaping data -------------------------------------------------------------
  nrow_origin <- nrow(x)

  # TODO
  # x, yにNAがある場合は、x,yを取り除く必要がある

  # Checking if x is a gridded data
  if (ncol(x) > 1) {

    is_gridded <- extract_gridded_x(x)

    x <- x[is_gridded, , drop = FALSE]
    y <- y[is_gridded]

    if (any(!is_gridded)) {
      cat(sum(!is_gridded), '/', nrow_origin,
          'rows of x were not gridded and dropped.\n')
    }

  }

  # Averaging y with duplicated x
  duplication_id <- detect_duplicated_x(x)
  id_uniq <- unique(duplication_id[!is.na(duplication_id)])
  if (length(id_uniq) > 0) {
    cat(sum(!is.na(duplication_id)), '/', nrow_origin,
        'were duplicated and their y values are averaged.\n')
  }
  dropped <- rep(FALSE, nrow(x))
  for (id in id_uniq) {

    id_index <- which(duplication_id == id)
    y[id_index[1]] <- mean(y[id_index])
    dropped[id_index[-1]] <- TRUE

  }
  x <- x[!dropped, , drop = FALSE]
  y <- y[!dropped]

  # TODO
  # 結果的にx,yがとても少なくなったらどうなるのか？
  # エラー文を出して処理を止める必要がある

  # Converting x to a list and y to an array
  x_list <- lapply(x, unique)
  y_arr <- tapply(X = y, INDEX = as.list(x), FUN = identity)
  names(dimnames(y_arr)) <- colnames(x)


  # Interpolation --------------------------------------------------------------

  # TODO
  # 日本語を英語に直す
  # rule = 1, 2の場合を実装する
  # verbose = TRUEの場合の記述内容を設定する

  # 入力された気象条件(dat)は、Nを入力されたデータ点数、
  # Kを気象条件の数(7)としたとき、N×Kの行列である。
  # dat内の要素をx_(n, k)と表記する。

  # x_(n, k)がrefdatの点のi番目とi+1番目にあるとき、そのiを計算
  # 答えはN×Kの行列(smaller_index)になる
  # if (verbose) t0 <- proc.time()
  # if (verbose) cat('\n  Fetching index... ')
  smaller_indices <- mapply(
    FUN = function(x, vals) {

      # iを計算
      # i <- sapply(vals, function(val) which(val - x < 0)[1] - 1)
      i <- sapply(vals, function(val) sum(val - x >= 0))

      # x_(n, k)がrefdatの点の範囲外にあるとき、
      # smaller_indexは0やlength(x)といった値を取るため、datこれを補正する。
      i[i == 0] <- 1
      i[i == length(x)] <- length(x) - 1
      return(i)

    },
    x = x_list,
    vals = xout)
  # if (verbose) cat('Done. (', round(proc.time() - t0, 2)[3], ' s)\n', sep = '')

  # x_(n, k)がrefdatの点のi番目(xr_(k, i))とi+1番目(xr_(k, i+1))にあるとき、
  # xr_(k, i+1) - xr_(k, i)を計算
  # 答えはN×Kの行列(region_length)になる
  # if (verbose) t0 <- proc.time()
  # if (verbose) cat('  Fetching region length... ')
  region_length <- mapply(
    FUN   = function(x, i) x[i + 1] - x[i],
    x = x_list,
    i     = as.data.frame(smaller_indices))
  # if (verbose) cat('Done. (', round(proc.time() - t0, 2)[3], ' s)\n', sep = '')

  # データ点を囲むrefdat内の格子点すべてに対して計算を行うため、
  # 格子点すべてのパターンを作っておく
  added_indices <- do.call(expand.grid, rep(list(0:1), ncol(x)))
  colnames(added_indices) <- colnames(x)

  # 答え(近似計算された穂温)を保存するベクトル(長さN)を用意しておく
  yout <- rep(0, nrow(xout))

  # データ点を囲む格子点一つずつに対し処理を行う
  # if (verbose) t0 <- proc.time()
  # if (verbose) cat('  Calculating weighted mean... ')
  for (j in seq_len(nrow(added_indices))) {
    # j <- 2

    added_index_j <- unlist(added_indices[j, ])

    # 現在注目している格子点の穂温を取得
    # t1 <- proc.time()
    y_indices_j <- smaller_indices + added_index_j[col(smaller_indices)]
    y_j <- y_arr[y_indices_j]
    # proc.time() - t1

    # 現在注目している格子点にかける重みを計算
    # t1 <- proc.time()
    rev_added_index_j <- 1 - added_index_j
    x_rev_indices_j <- smaller_indices + rev_added_index_j[col(smaller_indices)]
    sign_adjuster_j <- added_index_j * 2 - 1
    # proc.time() - t1
    # 0.08

    # t1 <- proc.time()
    weight_raw_j <- mapply(
      FUN = function(x, reglen, ref_x, x_rev_ind, sign_adj) {
        sign_adj * (x - ref_x[x_rev_ind]) / reglen
      },
      x         = xout,                            # N * 6
      reglen    = as.data.frame(region_length),   # N * 6
      ref_x     = x_list,                       # list(breaks * 6)
      x_rev_ind = as.data.frame(x_rev_indices_j), # N * 6
      sign_adj  = sign_adjuster_j)                # 6
    # proc.time() - t1
    # 0.50

    # t1 <- proc.time()
    # weight_j <- apply(X = weight_raw_j, MARGIN = 1, FUN = prod)
    # と同じことをやりたいのだが、apply関数は時間がかかるため、rowSumsで代替
    weight_j <- rep(0, nrow(weight_raw_j))
    is_0 <- rowSums(weight_raw_j == 0) > 0
    prodabs_weight_raw_j <- exp(rowSums(log(abs(weight_raw_j[!is_0, , drop = FALSE]))))
    n_minus_weight_raw_j <- rowSums(sign(weight_raw_j[!is_0, , drop = FALSE]) == -1)
    weight_j[!is_0] <- (-1) ^ (n_minus_weight_raw_j %% 2) * prodabs_weight_raw_j
    # proc.time() - t1
    # 0.58

    # 答えを更新
    yout <- yout + y_j * weight_j

  }
  # if (verbose) cat('Done. (', round(proc.time() - t0, 2)[3], ' s)\n', sep = '')
  return(yout)

}
