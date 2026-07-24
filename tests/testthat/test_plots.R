
set.seed(12345)
x <- sdetorus::toPiInt(mvtnorm::rmvnorm(n = 100, mean = c(0, 1),
                                        sigma = rbind(c(1, 0.3), c(0.3, 1))))

test_that("torus_pairs builds its diagonal and off-diagonal panels", {

  # Printing the ggmatrix forces evaluation of the panel closures
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  p <- torus_pairs(x)
  expect_s3_class(p, "ggmatrix")
  suppressWarnings(print(p))

  # Non-pi scales exercise the "l != pi" axis-label branches
  x2 <- sdetorus::toInt(x, a = -2, b = 2)
  suppressWarnings(print(torus_pairs(x2, scales = c(2, 2))))

  # Alternative bandwidth selectors exercise other switch() arms
  suppressWarnings(print(torus_pairs(x, bwd = "AMI")))
  suppressWarnings(print(torus_pairs(x, bwd = "EMI")))

})

test_that("show_ridge_pca covers the no-projection branch", {

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  fit <- ridge_pca(x = x, type = "bvm")
  expect_null(show_ridge_pca(fit, projs = FALSE, projs_lines = FALSE,
                             n_max = 40, signs = FALSE))
  expect_null(show_ridge_pca(fit, projs = TRUE, projs_lines = FALSE,
                             main = "ridge"))

})
