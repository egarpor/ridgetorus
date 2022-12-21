
set.seed(23131)
n <- 500
mu1 <- c(-pi, 1)
kappa1 <- c(2, 1, 0)
mu2 <- c(1, 2)
kappa2 <- c(3:1)
mu3 <- c(0, 0)
kappa3 <- c(3, 1, 5)
mu4 <- c(-1, -2)
kappa4 <- c(0.3, 1.1, 0.5)
mu5 <- c(-1, pi)
kappa5 <- c(0.1, 0.8, -2)
mu6 <- c(0, 1)
kappa6 <- c(0.5, 0.2, -2.5)
mu7 <- c(-3, -2)
kappa7 <- c(2, 2, -0.3)
mu8 <- c(0, -pi)
kappa8 <- c(1, 1, -5)
mu9 <- c(3, -2)
kappa9 <- c(1, 2, 0)
mu10 <- c(0, 2)
kappa10 <- c(10, 1, 0)
mu11 <- c(0, 2)
kappa11 <- c(5, 5, 0)
mu12 <- c(0, 2)
kappa12 <- c(10, 10, 0)

samp1 <- r_bvm(n = n, mu = mu1, kappa = kappa1)
samp2 <- r_bvm(n = n, mu = mu2, kappa = kappa2)
samp3 <- r_bvm(n = n, mu = mu3, kappa = kappa3)
samp4 <- r_bvm(n = n, mu = mu4, kappa = kappa4)
samp5 <- r_bvm(n = n, mu = mu5, kappa = kappa5)
samp6 <- r_bvm(n = n, mu = mu6, kappa = kappa6)
samp7 <- r_bvm(n = n, mu = mu7, kappa = kappa7)
samp8 <- r_bvm(n = n, mu = mu8, kappa = kappa8)
samp9 <- r_bvm(n = n, mu = mu9, kappa = kappa9)
samp10 <- r_bvm(n = n, mu = mu10, kappa = kappa10)
samp11 <- r_bvm(n = n, mu = mu11, kappa = kappa11)
samp12 <- r_bvm(n = n, mu = mu12, kappa = kappa12)

test_that("Integrates one in standard cases", {

  for (k1 in c(1, 5, 30)) {
    for (k2 in c(1, 5, 30)) {
      for (lamb in c(1, 5, 30)) {

        log_const <- log(const_bvm(kappa = c(k1, k2, lamb), M = 25))
        expect_equal(sdetorus::mcTorusIntegrate(
          f = d_bvm, mu = 1:2, kappa = c(k1, k2, lamb), p = 2,
          M = 1e5, log_const = log_const), 1, tolerance = 0.1)
        expect_equal(sdetorus::mcTorusIntegrate(
          f = d_bvm, mu = 1:2, kappa = c(k1, k2, lamb), p = 2,
          M = 1e5), 1, tolerance = 0.1)

      }
    }
  }

})

test_that("Integrates one for a concentration close to 0, for lambda <= 30", {

  for (lamb in c(0, 1, 15, 30)) {
    for (k2 in c(1, 15, 30)) {

      log_const <- suppressWarnings(
        log(const_bvm(kappa = c(1e-10, k2, lamb), M = 25)))
      expect_equal(sdetorus::mcTorusIntegrate(
        f = d_bvm, mu = 1:2, kappa = c(1e-10, k2, lamb),
        log_const = log_const, p = 2, M = 1e5), 1, tolerance = 0.1)

    }
  }

})

test_that("Integrates one for concentrations close to 0, for lambda <= 30", {

  for (lamb in c(0, 1, 15, 30)) {

    log_const <- suppressWarnings(
      log(const_bvm(kappa = c(1e-10, 1e-10, lamb))))
    expect_equal(sdetorus::mcTorusIntegrate(
      f = d_bvm, mu = 1:2, kappa = c(1e-10, 1e-10, lamb),
      log_const = log_const, p = 2, M = 1e5), 1, tolerance = 0.1)

  }

})

test_that("Integrates one for concentrations close to 0, for lambda > 30", {

  lamb <- 45
  log_const <- suppressWarnings(
    log(const_bvm(kappa = c(1e-10, 1e-10, lamb))))
  expect_equal(sdetorus::mcTorusIntegrate(
    f = d_bvm, mu = 1:2, kappa = c(1e-10, 1e-10, lamb),
    log_const = log_const, p = 2, M = 1e5), 1, tolerance = 0.1)

})

test_that("Recover parameters via MM", {

  fit1 <- unname(unlist(fit_bvm_mm(samp1)))
  fit2 <- unname(unlist(fit_bvm_mm(samp2)))
  fit3 <- unname(unlist(fit_bvm_mm(samp3)))
  fit4 <- unname(unlist(fit_bvm_mm(samp4)))
  fit5 <- unname(unlist(fit_bvm_mm(samp5)))
  fit6 <- unname(unlist(fit_bvm_mm(samp6, start = c(1, 0.5, 2, 0.5, 0.7))))

  expect_equal(c(torus_dist(x = rbind(fit1[1:2]), y = rbind(mu1[1:2])),
                 fit1[3:5] - kappa1), rep(0, 4), tolerance = 0.5)
  expect_equal(c(torus_dist(x = rbind(fit2[1:2]), y = rbind(mu2[1:2])),
                 fit2[3:5] - kappa2), rep(0, 4), tolerance = 0.5)
  expect_equal(c(torus_dist(x = rbind(fit3[1:2]), y = rbind(mu3[1:2])),
                 fit3[3:5] - kappa3), rep(0, 4), tolerance = 0.5)
  expect_equal(c(torus_dist(x = rbind(fit4[1:2]), y = rbind(mu4[1:2])),
                 fit4[3:5] - kappa4), rep(0, 4), tolerance = 0.5)
  expect_equal(c(torus_dist(x = rbind(fit5[1:2]), y = rbind(mu5[1:2])),
                 fit5[3:5] - kappa5), rep(0, 4), tolerance = 0.5)
  expect_equal(c(torus_dist(x = rbind(fit6[1:2]), y = rbind(mu6[1:2])),
                 fit6[3:5] - kappa6), rep(0, 4), tolerance = 0.5)

})

test_that("Recover parameters via MLE", {

  fit1 <- unname(unlist(fit_bvm_mle(samp1)[1:5]))
  fit2 <- unname(unlist(fit_bvm_mle(samp2)[1:5]))
  fit3 <- unname(unlist(fit_bvm_mle(samp3)[1:5]))
  fit4 <- unname(unlist(fit_bvm_mle(samp4)[1:5]))
  fit5 <- unname(unlist(fit_bvm_mle(samp5)[1:5]))
  fit6 <- unname(unlist(fit_bvm_mle(samp6)[1:5]))

  expect_equal(c(torus_dist(x = rbind(fit1[1:2]), y = rbind(mu1[1:2])),
                 fit1[3:5] - kappa1), rep(0, 4), tolerance = 0.5)
  expect_equal(c(torus_dist(x = rbind(fit2[1:2]), y = rbind(mu2[1:2])),
                 fit2[3:5] - kappa2), rep(0, 4), tolerance = 0.5)
  expect_equal(c(torus_dist(x = rbind(fit3[1:2]), y = rbind(mu3[1:2])),
                 fit3[3:5] - kappa3), rep(0, 4), tolerance = 0.5)
  expect_equal(c(torus_dist(x = rbind(fit4[1:2]), y = rbind(mu4[1:2])),
                 fit4[3:5] - kappa4), rep(0, 4), tolerance = 0.5)
  expect_equal(c(torus_dist(x = rbind(fit5[1:2]), y = rbind(mu5[1:2])),
                 fit5[3:5] - kappa5), rep(0, 4), tolerance = 0.5)
  expect_equal(c(torus_dist(x = rbind(fit6[1:2]), y = rbind(mu6[1:2])),
                 fit6[3:5] - kappa6), rep(0, 4), tolerance = 0.5)

})

test_that("Recover parameters under homogeneity via MM", {

  fit7 <- unname(unlist(fit_bvm_mm(samp7, hom = TRUE)))
  fit8 <- unname(unlist(fit_bvm_mm(samp8, hom = TRUE)))

  expect_equal(fit7[3], fit7[4])
  expect_equal(fit8[3], fit8[4])
  expect_equal(c(torus_dist(x = rbind(fit7[1:2]), y = rbind(mu7[1:2])),
                 fit7[3:5] - kappa7), rep(0, 4), tolerance = 0.5)
  expect_equal(c(torus_dist(x = rbind(fit8[1:2]), y = rbind(mu8[1:2])),
                 fit8[3:5] - kappa8), rep(0, 4), tolerance = 0.5)

})

test_that("Recover parameters under homogeneity via MLE", {

  fit7 <- unname(unlist(fit_bvm_mle(samp7, hom = TRUE)[1:5]))
  fit8 <- unname(unlist(fit_bvm_mle(samp8, hom = TRUE)[1:5]))

  expect_equal(fit7[3], fit7[4])
  expect_equal(fit8[3], fit8[4])
  expect_equal(c(torus_dist(x = rbind(fit7[1:2]), y = rbind(mu7[1:2])),
                 fit7[3:5] - kappa7), rep(0, 4), tolerance = 0.5)
  expect_equal(c(torus_dist(x = rbind(fit8[1:2]), y = rbind(mu8[1:2])),
                 fit8[3:5] - kappa8), rep(0, 4), tolerance = 0.5)

})

test_that("Recover parameters under independence via MM", {

  fit9 <- unname(unlist(fit_bvm_mm(samp9, indep = TRUE)[1:5]))
  fit10 <- unname(unlist(fit_bvm_mm(samp10, indep = TRUE)[1:5]))

  expect_equal(fit9[5], 0)
  expect_equal(fit10[5], 0)
  expect_equal(c(torus_dist(x = rbind(fit9[1:2]), y = rbind(mu9[1:2])),
                 fit9[3:5] - kappa9), rep(0, 4), tolerance = 0.5)
  expect_equal(c(torus_dist(x = rbind(fit10[1:2]), y = rbind(mu10[1:2])),
                 fit10[3:5] - kappa10), rep(0, 4), tolerance = 0.5)

})

test_that("Recover parameters under independence via MLE", {

  fit9 <- unname(unlist(fit_bvm_mle(samp9, indep = TRUE)[1:5]))
  fit10 <- unname(unlist(fit_bvm_mle(samp10, indep = TRUE)[1:5]))

  expect_equal(fit9[5], 0)
  expect_equal(fit10[5], 0)
  expect_equal(c(torus_dist(x = rbind(fit9[1:2]), y = rbind(mu9[1:2])),
                 fit9[3:5] - kappa9), rep(0, 4), tolerance = 0.5)
  expect_equal(c(torus_dist(x = rbind(fit10[1:2]), y = rbind(mu10[1:2])),
                 fit10[3:5] - kappa10), rep(0, 4), tolerance = 0.5)

})

test_that("Recover parameters under homogeneity and independence via MM", {

  fit11 <- unname(unlist(fit_bvm_mm(samp11, hom = TRUE, indep = TRUE)))
  fit12 <- unname(unlist(fit_bvm_mm(samp12, hom = TRUE, indep = TRUE)))

  expect_equal(fit11[5], 0)
  expect_equal(fit12[5], 0)
  expect_equal(c(torus_dist(x = rbind(fit11[1:2]), y = rbind(mu11[1:2])),
                 fit11[3:5] - kappa11), rep(0, 4), tolerance = 0.5)
  expect_equal(c(torus_dist(x = rbind(fit12[1:2]), y = rbind(mu12[1:2])),
                 fit12[3:5] - kappa12), rep(0, 4), tolerance = 0.5)

})

test_that("Recover parameters under homogeneity and independence via MLE", {

  fit11 <- unname(unlist(fit_bvm_mle(samp11, hom = TRUE, indep = TRUE)[1:5]))
  fit12 <- unname(unlist(fit_bvm_mle(samp12, hom = TRUE, indep = TRUE)[1:5]))

  expect_equal(fit11[5], 0)
  expect_equal(fit12[5], 0)
  expect_equal(c(torus_dist(x = rbind(fit11[1:2]), y = rbind(mu11[1:2])),
                 fit11[3:5] - kappa11), rep(0, 4), tolerance = 0.5)
  expect_equal(c(torus_dist(x = rbind(fit12[1:2]), y = rbind(mu12[1:2])),
                 fit12[3:5] - kappa12), rep(0, 4), tolerance = 0.5)

})

test_that("BvM errors", {

  expect_error(fit_bvm_mm(x = c(0, 0), start = 4))
  expect_error(fit_bvm_mle(x = c(0, 0), start = 4))

})
