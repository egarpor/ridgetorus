
set.seed(1)
# Sample from a wrapped normal
mu <- 2.3
std <- 1.4
x <- sdetorus::toPiInt(rnorm(n = 5e2, mean = mu, sd = std))

#Periodic data in [-2, 2)
mu2 <- 0
std2 <- 0.5
x2 <- sdetorus::toInt(rnorm(n = 5e2, mean = mu2, sd = std2), a = -2, b = 2)

test_that("frechet_mean estimates parameters correctly", {

  frechet <- frechet_mean(x = x, draw_plot = TRUE)
  frechet2 <- frechet_mean(x = x2, draw_plot = TRUE, l = 2)
  expect_equal(frechet$mu, mu, tolerance = 0.1)
  expect_equal(frechet$sd, std, tolerance = 0.1)
  expect_equal(frechet$var, std^2, tolerance = 0.2)
  expect_equal(frechet2$mu, mu2, tolerance = 0.1)
  expect_equal(frechet2$sd, std2, tolerance = 0.1)
  expect_equal(frechet2$var, std2^2, tolerance = 0.2)

})
