
# Common parameters
M_0 <- 200
n_0 <- 300
M_1 <- 50
n_1 <- 50
alpha <- c(0.05, 0.10, 0.25)
delta <- 0.75
mu <- c(0, 0)

# Homogeneity
kappa_01 <- c(1, 1, 0.25)
kappa_02 <- c(3, 3, -10)
kappa_11 <- c(1, 2, 0.25)
kappa_12 <- c(1, 5, 10)
xi_01 <- c(0.3, 0.3, 0.25)
xi_02 <- c(0.7, 0.7, 0.85)
xi_11 <- c(0.1, 0.5, 0.25)
xi_12 <- c(0.6, 0.2, 0.1)

# Independence
kappa_03 <- c(0.5, 1, 0)
kappa_04 <- c(3, 0.5, 0)
kappa_13 <- c(0.5, 1, 3)
kappa_14 <- c(2, 0.5, -1)
xi_03 <- c(0.2, 0.5, 0)
xi_04 <- c(0.8, 0.3, 0)
xi_13 <- c(0.2, 0.3, 0.2)
xi_14 <- c(0.8, 0.6, -0.2)

# Homogeneity and independence
kappa_05 <- c(5, 5, 0)
kappa_06 <- c(2, 2, 0)
kappa_15 <- c(3, 5, 0)
kappa_16 <- c(2, 2, 2)
xi_05 <- c(0.5, 0.5, 0)
xi_06 <- c(0.25, 0.25, 0)
xi_15 <- c(0, 0.3, 0.1)
xi_16 <- c(0.5, 0.1, 0.2)

## Homogeneity tests

set.seed(123456)

test_that("Homogeneity tests for bvm under H0", {

  p01 <- replicate(M_0, biv_lrt(x = r_bvm(n = n_0, mu = mu, kappa = kappa_01),
                                hom = TRUE, type = "bvm")$p.value)
  p02 <- replicate(M_0, biv_lrt(x = r_bvm(n = n_0, mu = mu, kappa = kappa_02),
                                hom = TRUE, type = "bvm")$p.value)
  expect_true(max(abs(unname(quantile(p01, probs = alpha)) / alpha - 1))
              < delta)
  expect_true(max(abs(unname(quantile(p02, probs = alpha)) / alpha - 1))
              < delta)

})

test_that("Homogeneity tests for bwc under H0", {

  p01 <- replicate(M_0, biv_lrt(x = r_bwc(n = n_0, mu = mu, xi = xi_01),
                                hom = TRUE, type = "bwc")$p.value)
  p02 <- replicate(M_0, biv_lrt(x = r_bwc(n = n_0, mu = mu, xi = xi_02),
                                hom = TRUE, type = "bwc")$p.value)
  expect_true(max(abs(unname(quantile(p01, probs = alpha)) / alpha - 1))
              < delta)
  expect_true(max(abs(unname(quantile(p02, probs = alpha)) / alpha - 1))
              < delta)

})

test_that("Homogeneity tests for bvm under H1", {

  p11 <- replicate(M_1, biv_lrt(x = r_bvm(n = n_1, mu = mu, kappa = kappa_11),
                                hom = TRUE, type = "bvm")$p.value)
  p12 <- replicate(M_1, biv_lrt(x = r_bvm(n = n_1, mu = mu, kappa = kappa_12),
                                hom = TRUE, type = "bvm")$p.value)
  expect_true(mean(p11) < 0.30)
  expect_true(mean(p12) < 0.30)

})

test_that("Homogeneity tests for bwc under H1", {

  p11 <- replicate(M_1, biv_lrt(x = r_bwc(n = n_1, mu = mu, xi = xi_11),
                               hom = TRUE, type = "bwc")$p.value)
  p12 <- replicate(M_1, biv_lrt(x = r_bwc(n = n_1, mu = mu, xi = xi_12),
                               hom = TRUE, type = "bwc")$p.value)
  expect_true(mean(p11) < 0.30)
  expect_true(mean(p12) < 0.30)

})

## Independence tests

set.seed(123456)

test_that("Independence tests for bvm under H0", {

  p01 <- replicate(M_0, biv_lrt(x = r_bvm(n = n_0, mu = mu, kappa = kappa_03),
                                indep = TRUE, type = "bvm")$p.value)
  p02 <- replicate(M_0, biv_lrt(x = r_bvm(n = n_0, mu = mu, kappa = kappa_04),
                                indep = TRUE, type = "bvm")$p.value)
  expect_true(max(abs(unname(quantile(p01, probs = alpha)) / alpha - 1))
              < delta)
  expect_true(max(abs(unname(quantile(p02, probs = alpha)) / alpha - 1))
              < delta)

})

test_that("Independence tests for bwc under H0", {

  p01 <- replicate(M_0, biv_lrt(x = r_bwc(n = n_0, mu = mu, xi = xi_03),
                                indep = TRUE, type = "bwc")$p.value)
  p02 <- replicate(M_0, biv_lrt(x = r_bwc(n = n_0, mu = mu, xi = xi_04),
                                indep = TRUE, type = "bwc")$p.value)
  expect_true(max(abs(unname(quantile(p01, probs = alpha)) / alpha - 1))
              < delta)
  expect_true(max(abs(unname(quantile(p02, probs = alpha)) / alpha - 1))
              < delta)

})

test_that("Independence tests for bvm under H1", {

  p11 <- replicate(M_1, biv_lrt(x = r_bvm(n = n_1, mu = mu, kappa = kappa_13),
                                indep = TRUE, type = "bvm")$p.value)
  p12 <- replicate(M_1, biv_lrt(x = r_bvm(n = n_1, mu = mu, kappa = kappa_14),
                                indep = TRUE, type = "bvm")$p.value)
  expect_true(mean(p11) < 0.30)
  expect_true(mean(p12) < 0.30)

})

test_that("Independence tests for bwc under H1", {

  p11 <- replicate(M_1, biv_lrt(x = r_bwc(n = n_1, mu = mu, xi = xi_13),
                                indep = TRUE, type = "bwc")$p.value)
  p12 <- replicate(M_1, biv_lrt(x = r_bwc(n = n_1, mu = mu, xi = xi_14),
                                indep = TRUE, type = "bwc")$p.value)
  expect_true(mean(p11) < 0.30)
  expect_true(mean(p12) < 0.30)

})

## Homogeneity and independence tests

set.seed(123456)

test_that("Homogeneity and independence for bvm under H0", {

  p01 <- replicate(M_0, biv_lrt(x = r_bvm(n = n_0, mu = mu, kappa = kappa_05),
                                indep = TRUE, hom = TRUE, type = "bvm")$p.value)
  p02 <- replicate(M_0, biv_lrt(x = r_bvm(n = n_0, mu = mu, kappa = kappa_06),
                                indep = TRUE, hom = TRUE, type = "bvm")$p.value)
  expect_true(max(abs(unname(quantile(p01, probs = alpha)) / alpha - 1))
              < delta)
  expect_true(max(abs(unname(quantile(p02, probs = alpha)) / alpha - 1))
              < delta)

})

test_that("Homogeneity and independence for bwc under H0", {

  p01 <- replicate(M_0, biv_lrt(x = r_bwc(n = n_0, mu = mu, xi = xi_05),
                                indep = TRUE, hom = TRUE, type = "bwc")$p.value)
  p02 <- replicate(M_0, biv_lrt(x = r_bwc(n = n_0, mu = mu, xi = xi_06),
                                indep = TRUE, hom = TRUE, type = "bwc")$p.value)
  expect_true(max(abs(unname(quantile(p01, probs = alpha)) / alpha - 1))
              < delta)
  expect_true(max(abs(unname(quantile(p02, probs = alpha)) / alpha - 1))
              < delta)

})

test_that("Homogeneity and independence for bvm under H1", {

  p11 <- replicate(M_1, biv_lrt(x = r_bvm(n = n_1, mu = mu, kappa = kappa_15),
                                indep = TRUE, hom = TRUE, type = "bvm")$p.value)
  p12 <- replicate(M_1, biv_lrt(x = r_bvm(n = n_1, mu = mu, kappa = kappa_16),
                                indep = TRUE, hom = TRUE, type = "bvm")$p.value)
  expect_true(mean(p11) < 0.30)
  expect_true(mean(p12) < 0.30)

})

test_that("LRT error", {

  expect_error(biv_lrt(x = r_bwc(n = n_1, mu = mu, xi = xi_15), hom = FALSE,
                       indep = FALSE))

})

test_that("Homogeneity and independence for bwc under H1", {

  p11 <- replicate(M_1, biv_lrt(x = r_bwc(n = n_1, mu = mu, xi = xi_15),
                                indep = TRUE, hom = TRUE, type = "bwc")$p.value)
  p12 <- replicate(M_1, biv_lrt(x = r_bwc(n = n_1, mu = mu, xi = xi_16),
                                indep = TRUE, hom = TRUE, type = "bwc")$p.value)
  expect_true(mean(p11) < 0.30)
  expect_true(mean(p12) < 0.30)

})
