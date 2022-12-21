
# Samples
set.seed(12314)
n <- 200
samp_bvm_1 <- r_bvm(n = n, mu = 3:2, kappa = c(3, 1, 2))
samp_bvm_1_NA <- rbind(samp_bvm_1, c(NA, NA))
samp_bvm_2 <- r_bvm(n = n, mu = 3:2, kappa = c(2, 2, 5))
samp_bvm_3 <- r_bvm(n = n, mu = 3:2, kappa = c(5, 1, 0))
samp_bvm_4 <- r_bvm(n = n, mu = 3:2, kappa = c(1, 2, 0))
samp_bvm_5 <- r_bvm(n = n, mu = 3:2, kappa = c(2, 2, 0))
samp_bwc_1 <- r_bwc(n = n, mu = 3:2, xi = c(0.2, 0.8, 0.3))
samp_bwc_2 <- r_bwc(n = n, mu = 3:2, xi = c(0.2, 0.2, 0.5))
samp_bwc_3 <- r_bwc(n = n, mu = 3:2, xi = c(0.5, 0.1, 0))
samp_bwc_4 <- r_bwc(n = n, mu = 3:2, xi = c(0.2, 0.6, 0))
samp_bwc_5 <- r_bwc(n = n, mu = 3:2, xi = c(0.2, 0.2, 0))
samp_clust <- sdetorus::toPiInt(rbind(
                         mvtnorm::rmvnorm(n = n, mean = c(-pi, -pi) / 2,
                                          sigma = diag(0.1, nrow = 2)),
                         mvtnorm::rmvnorm(n = n, mean = c(-3 * pi / 2, 0) / 2,
                                 sigma = diag(0.1, nrow = 2)),
                         mvtnorm::rmvnorm(n = n, mean = c(0, pi / 2),
                                          sigma = diag(0.1, nrow = 2))))

# Ridges
fit_bvm_1 <- ridge_pca(x = samp_bvm_1)
fit_bvm_1_NA <- ridge_pca(x = samp_bvm_1_NA)
fit_bvm_2 <- ridge_pca(x = samp_bvm_2)
fit_bvm_3 <- ridge_pca(x = samp_bvm_3)
fit_bvm_4 <- ridge_pca(x = samp_bvm_4)
fit_bwc_1 <- ridge_pca(x = samp_bwc_1)
fit_bwc_2 <- ridge_pca(x = samp_bwc_2)
fit_bwc_3 <- ridge_pca(x = samp_bwc_3)
fit_bwc_4 <- ridge_pca(x = samp_bwc_4)
fit_clust <- ridge_pca(x = samp_clust, type = "bvm")

test_that("TR-PCA identifies correct model", {

  expect_equal(fit_bvm_1$type, "bvm")
  expect_equal(fit_bvm_1_NA$type, "bvm")
  expect_equal(fit_bvm_2$type, "bvm")
  expect_equal(fit_bvm_3$type, "bvm")
  expect_equal(fit_bvm_4$type, "bvm")
  expect_equal(fit_bwc_1$type, "bwc")
  expect_equal(fit_bwc_2$type, "bwc")
  expect_equal(fit_bwc_3$type, "bwc")
  expect_equal(fit_bwc_4$type, "bwc")

})

test_that("TR-PCA identifies regular/diagonal/horizontal/vertical ridges", {

  # Regular
  expect_true(all(fit_bvm_1$fit_mle[1:5] != 0))
  expect_true(all(fit_bwc_1$fit_mle[1:5] != 0))

  # Diagonal
  expect_equal(fit_bvm_2$fit_mle$kappa1, fit_bvm_2$fit_mle$kappa2)
  expect_equal(fit_bwc_2$fit_mle$xi1, fit_bwc_2$fit_mle$xi2)

  # Horizontal
  expect_equal(fit_bvm_3$fit_mle$lambda, 0)
  expect_gt(fit_bvm_3$fit_mle$kappa1, fit_bvm_3$fit_mle$kappa2)
  expect_equal(fit_bwc_3$fit_mle$rho, 0)
  expect_gt(fit_bwc_3$fit_mle$xi1, fit_bwc_3$fit_mle$xi2)

  # Vertical
  expect_equal(fit_bvm_4$fit_mle$lambda, 0)
  expect_gt(fit_bvm_4$fit_mle$kappa2, fit_bvm_4$fit_mle$kappa1)
  expect_equal(fit_bwc_4$fit_mle$rho, 0)
  expect_gt(fit_bwc_4$fit_mle$xi2, fit_bwc_4$fit_mle$xi1)

})

test_that("TR-PCA detects bvm problematic ridges", {

  expect_warning(fit_bvm_5 <- ridge_pca(x = samp_bvm_5))
  expect_true(all(fit_bvm_5$fit_mle[1:5] != 0))
  torus_pairs(x = samp_bvm_5)
  show_ridge_pca(fit_bvm_5, col_data = c(1, 2), signs = FALSE)

})

test_that("TR-PCA detects bwc problematic ridges", {

  expect_warning(fit_bwc_5 <- ridge_pca(x = samp_bwc_5))
  expect_true(all(fit_bwc_5$fit_mle[1:5] != 0))
  torus_pairs(x = samp_bwc_5)
  show_ridge_pca(fit_bwc_5)

})

test_that("ridge_pca errors", {

  expect_error(ridge_pca(x = c(1, 2, 3)))
  expect_error(ridge_pca(x = cbind(samp_bvm_1, samp_bvm_1)))
  expect_error(ridge_pca(x = samp_bvm_4, type = "other density"))
  expect_error(ridge_pca(x = samp_bvm_4, N = -1))
  expect_error(ridge_pca(x = samp_bvm_4, K = -1))
  expect_error(ridge_pca(x = samp_bvm_4, alpha = -1))
  expect_error(ridge_pca(x = samp_bvm_4, alpha = 10))
  expect_error(ridge_pca(x = samp_bvm_4, scale = 3))
  expect_error(ridge_pca(x = samp_bvm_4, lrts = 0))
  expect_error(ridge_pca(x = samp_bvm_4, at2 = 0))

})
