library(mapprox)

arith_seq <- function(a12, n) a12[1] + (a12[2] - a12[1]) * (n - 1)

expect_equal_randcpp <- function(..., element = NULL, expected) {

  if (is.null(element)) {
    expect_equal(mapprox(..., use_cpp = FALSE), expected)
    expect_equal(mapprox(..., use_cpp = TRUE), expected)
  } else {
    expect_equal(mapprox(..., use_cpp = FALSE)[[element]], expected)
    expect_equal(mapprox(..., use_cpp = TRUE) [[element]], expected)
  }

}

test_that('正常な入力で動作するか: 1変量', {

  x <- 1:10
  y <- c(1:3, 8:10, 7:4)
  xout <- seq(0, 11, 0.5)

  res_approx_1 <- approx(x, y, xout, rule = 1)$y
  res_approx_2 <- approx(x, y, xout, rule = 2)$y
  res_approx_3 <- res_approx_1
  res_approx_3[2:1]   <- arith_seq(res_approx_3[4:3], 3:4)
  res_approx_3[22:23] <- arith_seq(res_approx_3[20:21], 3:4)

  expect_equal_randcpp(x, y, xout, rule = 1, element = 'x', expected = as.data.frame(as.matrix(x)))
  expect_equal_randcpp(x, y, xout, rule = 1, element = 'y', expected = y)

  expect_equal_randcpp(x, y, xout, rule = 1, element = 'yout', expected = res_approx_1)
  expect_equal_randcpp(x, y, xout, rule = 2, element = 'yout', expected = res_approx_2)
  expect_equal_randcpp(x, y, xout, rule = 3, element = 'yout', expected = res_approx_3)

})
test_that('正常な入力で動作するか: 2変量', {

  f <- function(V1, V2) (V1 - 11.5) ^ 2 + V2 * 5
  x_v1 <- 11:15
  x <- expand.grid(V1 = x_v1, V2 = 101:103)
  y <- with(x, f(V1, V2))

  xout_v1 <- seq(11, 15, 0.5)
  xout_v2 <- c(101.5, 102.8, 104)
  xout <- expand.grid(V1 = xout_v1, V2 = xout_v2)

  yout <- by(xout, xout$V2, function(df) {
    approx(x_v1, f(x_v1, unique(df$V2)), df$V1)$y
  })
  yout <- matrix(rep(unlist(yout), 3), ncol = 3)
  sel <- xout$V2 == 104
  yout[sel, 1] <- NA
  yout[sel, 2] <- approx(x_v1, f(x_v1, 103), xout$V1[sel])$y

  expect_equal_randcpp(x, y, xout, rule = 1, element = 'x', expected = x)
  expect_equal_randcpp(x, y, xout, rule = 1, element = 'y', expected = y)

  expect_equal_randcpp(x, y, xout, rule = 1, element = 'yout', expected = yout[, 1])
  expect_equal_randcpp(x, y, xout, rule = 2, element = 'yout', expected = yout[, 2])
  expect_equal_randcpp(x, y, xout, rule = 3, element = 'yout', expected = yout[, 3])

  sampleid <- sample(seq(nrow(x)))
  x_shuf <- x[sampleid, ]
  y_shuf <- y[sampleid]
  expect_equal(mapprox(x, y, xout, rule = 3)$yout,
               mapprox(x_shuf, y_shuf, xout, rule = 3)$yout)
  expect_equal(mapprox(x, y, xout, rule = 3, use_cpp = TRUE)$yout,
               mapprox(x_shuf, y_shuf, xout, rule = 3, use_cpp = TRUE)$yout)

})
test_that('NAや格子状にないデータを弾けるか', {

  # correct data format
  x <- expand.grid(V1 = 11:15, V2 = 101:103)
  y <- with(x, (V1 - 12) ^ 2 + V2 * 5)
  xout <- expand.grid(V1 = seq(11, 16, 0.1), V2 = 101:103)
  res_correct <- mapprox(x, y, xout, rule = 1)

  # NA in x
  x1 <- rbind(x, data.frame(V1 = NA, V2 = 103))
  y1 <- c(y, 1000)
  expect_warning(mapprox(x1, y1, xout, rule = 1),
                 '1 NAs in x and 0 NAs in y')
  expect_equal(suppressWarnings(mapprox(x1, y1, xout, rule = 1)),
               res_correct)

  # non-gridded value
  x2 <- rbind(x, data.frame(V1 = 16, V2 = 103))
  y2 <- c(y, 1000)
  expect_warning(mapprox(x2, y2, xout, rule = 1),
                 '1 values of x were not gridded')
  expect_equal(suppressWarnings(mapprox(x2, y2, xout, rule = 1)),
               res_correct)

  # NA in y
  x3 <- rbind(x, data.frame(V1 = 11, V2 = 104))
  y3 <- c(y, NA)
  expect_warning(mapprox(x3, y3, xout, rule = 1),
                 '0 NAs in x and 1 NAs in y')
  expect_equal(suppressWarnings(mapprox(x3, y3, xout, rule = 1)),
               res_correct)

  # NA in xout
  xout4 <- rbind(xout, data.frame(V1 = NA, V2 = 104))
  expect_equal(mapprox(x, y, xout4, rule = 1)$yout, c(res_correct$yout, NA))

})
test_that('xが重複している場合、対応するyを平均化できるか', {

  # correct data format
  f <- function(V1, V2) (V1 - 12) ^ 2 + V2 * 5
  x <- expand.grid(V1 = 11:15, V2 = 101:103)
  y <- with(x, f(V1, V2))
  xout <- expand.grid(V1 = seq(11, 16, 0.5), V2 = 101:103)

  # duplicated value
  x <- rbind(x, data.frame(V1 = 14, V2 = 103))
  y <- c(y, 1000)

  # correct answer
  yout <- NULL
  for (v2 in 101:103) {
    y_i <- f(11:15, v2)
    if (v2 == 103) y_i[4] <- mean(c(f(14, v2), 1000))
    yout_i <- approx(11:15, y_i, xout = seq(11, 16, 0.5))$y
    yout <- c(yout, yout_i)
  }

  # test
  expect_warning(mapprox(x, y, xout, rule = 1),
                 '2 values of x were duplicated')

})
test_that('不均一な格子点のデータ', {

  f <- function(V1, V2) (V1 - 11.5) ^ 2 + V2 * 5
  x_v1 <- c(11, 12, 15, 16)
  x <- expand.grid(V1 = x_v1, V2 = c(101, 101.5, 104))
  y <- with(x, f(V1, V2))

  xout_v1 <- seq(11, 16, 0.5)
  xout_v2 <- c(101.2, 103, 104)
  xout <- expand.grid(V1 = xout_v1, V2 = xout_v2)

  yout <- by(xout, xout$V2, function(df) {
    approx(x_v1, f(x_v1, unique(df$V2)), df$V1, rule = 2)$y
  })
  yout <- matrix(rep(unlist(yout), 3), ncol = 3)
  sel <- xout$V2 == 104.5
  yout[sel, 1] <- NA
  yout[sel, 2] <- approx(x_v1, f(x_v1, 103), xout$V1[sel])$y

  expect_equal(mapprox(x, y, xout, rule = 1)$x, x)
  expect_equal(mapprox(x, y, xout, rule = 1)$y, y)

  expect_equal(mapprox(x, y, xout, rule = 1)$yout, yout[, 1])
  expect_equal(mapprox(x, y, xout, rule = 2)$yout, yout[, 2])
  expect_equal(mapprox(x, y, xout, rule = 3)$yout, yout[, 3])

})
test_that('使用できるデータが極端に少ない場合にどうなるか', {

  # 使用できるx, yが０の場合
  # 使用できるxoutが０の場合

  # xのある変数において、使用できる値が１の場合
  # -> 未実装、警告を出して変数を消すのがよいだろう
  # -> その場合、１変数の場合は？
  # -> xが１つしかないということなので、エラーでよさそう

})
