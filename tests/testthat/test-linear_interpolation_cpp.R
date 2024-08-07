library(mapprox)

arith_seq <- function(a12, n) a12[1] + (a12[2] - a12[1]) * (n - 1)

test_that('正常な入力で動作するか: 1変量', {

  x <- data.frame(1:10)
  y <- c(1:3, 8:10, 7:4)
  xout <- data.frame(seq(0, 11, 0.5))

  yout <- linear_interpolation_cpp(x, y, xout, rule = 2, verbose = TRUE)
  yout

  # res_approx_1 <- approx(x, y, xout, rule = 1)$y
  # res_approx_2 <- approx(x, y, xout, rule = 2)$y
  # res_approx_3 <- res_approx_1
  # res_approx_3[2:1]   <- arith_seq(res_approx_3[4:3], 3:4)
  # res_approx_3[22:23] <- arith_seq(res_approx_3[20:21], 3:4)
  #
  # expect_equal(mapprox(x, y, xout, rule = 1)$x, as.data.frame(as.matrix(x)))
  # expect_equal(mapprox(x, y, xout, rule = 1)$y, y)
  #
  # expect_equal(mapprox(x, y, xout, rule = 1)$yout, res_approx_1)
  # expect_equal(mapprox(x, y, xout, rule = 2)$yout, res_approx_2)
  # expect_equal(mapprox(x, y, xout, rule = 3)$yout, res_approx_3)

})
test_that('正常な入力で動作するか: 2変量', {

  f <- function(V1, V2) (V1 - 11.5) ^ 2 + V2 * 5
  x_v1 <- 11:15
  x <- expand.grid(V1 = rev(x_v1), V2 = 101:103)
  y <- with(x, f(V1, V2))

  xout_v1 <- seq(10, 16, 0.5)
  xout_v2 <- c(101.5, 102.8, 104)
  xout <- expand.grid(V1 = xout_v1, V2 = xout_v2)

  yout <- linear_interpolation_cpp(x, y, xout, rule = 2, verbose = TRUE)
  yout

  # yout <- by(xout, xout$V2, function(df) {
  #   approx(x_v1, f(x_v1, unique(df$V2)), df$V1)$y
  # })
  # yout <- matrix(rep(unlist(yout), 3), ncol = 3)
  # sel <- xout$V2 == 104
  # yout[sel, 1] <- NA
  # yout[sel, 2] <- approx(x_v1, f(x_v1, 103), xout$V1[sel])$y
  #
  # expect_equal(mapprox(x, y, xout, rule = 1)$x, x)
  # expect_equal(mapprox(x, y, xout, rule = 1)$y, y)
  #
  # expect_equal(mapprox(x, y, xout, rule = 1)$yout, yout[, 1])
  # expect_equal(mapprox(x, y, xout, rule = 2)$yout, yout[, 2])
  # expect_equal(mapprox(x, y, xout, rule = 3)$yout, yout[, 3])

})
