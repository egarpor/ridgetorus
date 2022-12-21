
set.seed(12312)
n <- 500
mu1 <- c(-pi, 1)
Sigma1 <- matrix(c(1, 0, 0, 1), nrow = 2)
mu2 <- c(1, 2)
Sigma2 <- matrix(c(2, 0.5, 0.5, 1), nrow = 2)
mu3 <- c(0, 0)
Sigma3 <- matrix(c(3, 1, 1, 1), nrow = 2)
samp1 <- r_bwn(n = n, mu = mu1, Sigma = Sigma1)
samp2 <- r_bwn(n = n, mu = mu2, Sigma = Sigma2)
samp3 <- r_bwn(n = n, mu = mu3, Sigma = Sigma3)

test_that("Integrates one", {

  expect_equal(
    sdetorus::mcTorusIntegrate(f = d_bwn, mu = mu1, Sigma = Sigma1,
                               p = 2, M = 1e5), 1, tolerance = 0.01)
  expect_equal(
    sdetorus::mcTorusIntegrate(f = d_bwn, mu = mu1, Sigma = Sigma2,
                               p = 2, M = 1e5), 1, tolerance = 0.01)

})

test_that("Recover parameters via MLE", {

  fit1 <- fit_bwn_mle(samp1)
  fit2 <- fit_bwn_mle(samp2)
  fit3 <- fit_bwn_mle(samp3)# Cheating

  expect_equal(c(unname(torus_dist(x = rbind(fit1$mu[1]), y = mu1[1])),
                 unname(torus_dist(x = rbind(fit1$mu[2]), y = mu1[2])),
                 unname(fit1$Sigma[1, 1]) - Sigma1[1, 1],
                 unname(fit1$Sigma[1, 2]) - Sigma1[1, 2],
                 unname(fit1$Sigma[2, 1]) - Sigma1[2, 1],
                 unname(fit1$Sigma[2, 2]) - Sigma1[2, 2]),
               rep(0, 6), tolerance = 0.2)
  expect_equal(c(unname(torus_dist(x = rbind(fit2$mu[1]), y = mu2[1])),
                 unname(torus_dist(x = rbind(fit2$mu[2]), y = mu2[2])),
                 unname(fit2$Sigma[1, 1]) - Sigma2[1, 1],
                 unname(fit2$Sigma[1, 2]) - Sigma2[1, 2],
                 unname(fit2$Sigma[2, 1]) - Sigma2[2, 1],
                 unname(fit2$Sigma[2, 2]) - Sigma2[2, 2]),
               rep(0, 6), tolerance = 0.2)
  expect_equal(c(unname(torus_dist(x = rbind(fit3$mu[1]), y = mu3[1])),
                 unname(torus_dist(x = rbind(fit3$mu[2]), y = mu3[2])),
                 unname(fit3$Sigma[1, 1]) - Sigma3[1, 1],
                 unname(fit3$Sigma[1, 2]) - Sigma3[1, 2],
                 unname(fit3$Sigma[2, 1]) - Sigma3[2, 1],
                 unname(fit3$Sigma[2, 2]) - Sigma3[2, 2]),
               rep(0, 6), tolerance = 0.2)

})
