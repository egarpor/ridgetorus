
ridge <- ridge_bvm(mu = c(0, 0), kappa = c(0.3, 0.9, 0.8), subint_1 = 1e3,
                   subint_2 = 1e3)
coef <- ridge_fourier_fit(ridge, at2 = FALSE)
R_star <- function(phi, coef) {
  s <- 0
  for (i in seq_along(coef)) {
    s <- s + coef[i] * sin(i * phi)
  }
  return(s)
}
val <- R_star(phi = ridge[, 1], coef = coef)

test_that("The fit reproduces the ridge", {

  expect_equal(max(ridge[, 2] - val), 0, tolerance = 0.1)

})
