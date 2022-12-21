
set.seed(42341)
n <- 500
mu1 <- c(-pi, 1)
xi1 <- c(0.4, 0.7, 0)
mu2 <- c(1, 2)
xi2 <- c(0.5, 0.2, 0.4)
mu3 <- c(0, 0)
xi3 <- c(0.1, 0.6, 0.8)
mu4 <- c(-2, -1)
xi4 <- c(0.4, 0.7, -0.3)
mu5 <- c(2, 2.8)
xi5 <- c(0.5, 0.1, -0.7)
mu6 <- c(0, 0)
xi6 <- c(0.4, 0.3, 0)
mu7 <- c(2, 2.8)
xi7 <- c(0.25, 0.25, -0.9)
mu8 <- c(0, 0)
xi8 <- c(0.9, 0.9, 0.1)
mu9 <- c(1, 2)
xi9 <- c(0.1, 0.5, 0)
mu10 <- c(1, 2)
xi10 <- c(0.3, 0.5, 0)
mu11 <- c(1.5, 2)
xi11 <- c(0.7, 0.7, 0)
mu12 <- c(2, 1.5)
xi12 <- c(0.3, 0.3, 0)

samp1 <- r_bwc(n = n, mu = mu1, xi = xi1)
samp2 <- r_bwc(n = n, mu = mu2, xi = xi2)
samp3 <- r_bwc(n = n, mu = mu3, xi = xi3)
samp4 <- r_bwc(n = n, mu = mu4, xi = xi4)
samp5 <- r_bwc(n = n, mu = mu5, xi = xi5)
samp6 <- r_bwc(n = n, mu = mu6, xi = xi6)
samp7 <- r_bwc(n = n, mu = mu7, xi = xi7)
samp8 <- r_bwc(n = n, mu = mu8, xi = xi8)
samp9 <- r_bwc(n = n, mu = mu9, xi = xi9)
samp10 <- r_bwc(n = n, mu = mu10, xi = xi10)
samp11 <- r_bwc(n = n, mu = mu11, xi = xi11)
samp12 <- r_bwc(n = n, mu = mu12, xi = xi12)

test_that("Integrates one", {

  for (xii1 in seq(0, 0.8, l = 3)) {
    for (xii2 in seq(0, 0.8, l = 3)) {
      for (rho in seq(0, 0.8, l = 3)) {

      expect_equal(sdetorus::mcTorusIntegrate(
        f = d_bwc, mu = c(0, 0), xi = c(xii1, xii2, rho), p = 2, M = 1e5),
        1, tolerance = 0.1)

      }
    }
  }

})

test_that("Recover parameters via MM", {

  fit1 <- unname(unlist(fit_bwc_mm(samp1)))
  fit2 <- unname(unlist(fit_bwc_mm(samp2)))
  fit3 <- unname(unlist(fit_bwc_mm(samp3)))
  fit4 <- unname(unlist(fit_bwc_mm(samp4)))
  fit5 <- unname(unlist(fit_bwc_mm(samp5)))
  fit6 <- unname(unlist(fit_bwc_mm(samp6)))

  expect_equal(c(torus_dist(x = rbind(fit1[1:2]), y = rbind(mu1[1:2])),
                 fit1[3:5] - xi1), rep(0, 4), tolerance = 0.5)
  expect_equal(c(torus_dist(x = rbind(fit2[1:2]), y = rbind(mu2[1:2])),
                 fit2[3:5] - xi2), rep(0, 4), tolerance = 0.5)
  expect_equal(c(torus_dist(x = rbind(fit3[1:2]), y = rbind(mu3[1:2])),
                 fit3[3:5] - xi3), rep(0, 4), tolerance = 0.5)
  expect_equal(c(torus_dist(x = rbind(fit4[1:2]), y = rbind(mu4[1:2])),
                 fit4[3:5] - xi4), rep(0, 4), tolerance = 0.5)
  expect_equal(c(torus_dist(x = rbind(fit5[1:2]), y = rbind(mu5[1:2])),
                 fit5[3:5] - xi5), rep(0, 4), tolerance = 0.5)
  expect_equal(c(torus_dist(x = rbind(fit6[1:2]), y = rbind(mu6[1:2])),
                 fit6[3:5] - xi6), rep(0, 4), tolerance = 0.5)

})

test_that("Recover parameters via MLE", {

  fit1 <- unname(unlist(fit_bwc_mle(samp1)[1:5]))
  fit2 <- unname(unlist(fit_bwc_mle(samp2)[1:5]))
  fit3 <- unname(unlist(fit_bwc_mle(samp3)[1:5]))
  fit4 <- unname(unlist(fit_bwc_mle(samp4)[1:5]))
  fit5 <- unname(unlist(fit_bwc_mle(samp5)[1:5]))
  fit6 <- unname(unlist(fit_bwc_mle(samp6)[1:5]))

  expect_equal(c(torus_dist(x = rbind(fit1[1:2]), y = rbind(mu1[1:2])),
                 fit1[3:5] - xi1), rep(0, 4), tolerance = 0.5)
  expect_equal(c(torus_dist(x = rbind(fit2[1:2]), y = rbind(mu2[1:2])),
                 fit2[3:5] - xi2), rep(0, 4), tolerance = 0.5)
  expect_equal(c(torus_dist(x = rbind(fit3[1:2]), y = rbind(mu3[1:2])),
                 fit3[3:5] - xi3), rep(0, 4), tolerance = 0.5)
  expect_equal(c(torus_dist(x = rbind(fit4[1:2]), y = rbind(mu4[1:2])),
                 fit4[3:5] - xi4), rep(0, 4), tolerance = 0.5)
  expect_equal(c(torus_dist(x = rbind(fit5[1:2]), y = rbind(mu5[1:2])),
                 fit5[3:5] - xi5), rep(0, 4), tolerance = 0.5)
  expect_equal(c(torus_dist(x = rbind(fit6[1:2]), y = rbind(mu6[1:2])),
                 fit6[3:5] - xi6), rep(0, 4), tolerance = 0.5)

})

test_that("Recover parameters under homogeneity via MM", {

  fit7 <- unname(unlist(fit_bwc_mm(samp7, hom = TRUE)))
  fit8 <- unname(unlist(fit_bwc_mm(samp8, hom = TRUE)))

  expect_equal(fit7[3], fit7[4])
  expect_equal(fit8[3], fit8[4])
  expect_equal(c(torus_dist(x = rbind(fit7[1:2]), y = rbind(mu7[1:2])),
                 fit7[3:5] - xi7), rep(0, 4), tolerance = 0.5)
  expect_equal(c(torus_dist(x = rbind(fit8[1:2]), y = rbind(mu8[1:2])),
                 fit8[3:5] - xi8), rep(0, 4), tolerance = 0.5)

})

test_that("Recover parameters under homogeneity via MLE", {

  fit7 <- unname(unlist(fit_bwc_mle(samp7, hom = TRUE)[1:5]))
  fit8 <- unname(unlist(fit_bwc_mle(samp8, hom = TRUE)[1:5]))

  expect_equal(fit7[3], fit7[4])
  expect_equal(fit8[3], fit8[4])
  expect_equal(c(torus_dist(x = rbind(fit7[1:2]), y = rbind(mu7[1:2])),
                 fit7[3:5] - xi7), rep(0, 4), tolerance = 0.5)
  expect_equal(c(torus_dist(x = rbind(fit8[1:2]), y = rbind(mu8[1:2])),
                 fit8[3:5] - xi8), rep(0, 4), tolerance = 0.5)

})

test_that("Recover parameters under independence via MM", {

  fit9 <- unname(unlist(fit_bwc_mm(samp9, indep = TRUE)[1:5]))
  fit10 <- unname(unlist(fit_bwc_mm(samp10, indep = TRUE)[1:5]))

  expect_equal(fit9[5], 0)
  expect_equal(fit10[5], 0)
  expect_equal(c(torus_dist(x = rbind(fit9[1:2]), y = rbind(mu9[1:2])),
                 fit9[3:5] - xi9), rep(0, 4), tolerance = 0.5)
  expect_equal(c(torus_dist(x = rbind(fit10[1:2]), y = rbind(mu10[1:2])),
                 fit10[3:5] - xi10), rep(0, 4), tolerance = 0.5)

})

test_that("Recover parameters under independence via MLE", {

  fit9 <- unname(unlist(fit_bwc_mle(samp9, indep = TRUE)[1:5]))
  fit10 <- unname(unlist(fit_bwc_mle(samp10, indep = TRUE)[1:5]))

  expect_equal(fit9[5], 0)
  expect_equal(fit10[5], 0)
  expect_equal(c(torus_dist(x = rbind(fit9[1:2]), y = rbind(mu9[1:2])),
                 fit9[3:5] - xi9), rep(0, 4), tolerance = 0.5)
  expect_equal(c(torus_dist(x = rbind(fit10[1:2]), y = rbind(mu10[1:2])),
                 fit10[3:5] - xi10), rep(0, 4), tolerance = 0.5)

})

test_that("Recover parameters under homogeneity and independence via MM", {

  fit11 <- unname(unlist(fit_bwc_mm(samp11, hom = TRUE, indep = TRUE)))
  fit12 <- unname(unlist(fit_bwc_mm(samp12, hom = TRUE, indep = TRUE)))

  expect_equal(fit11[5], 0)
  expect_equal(fit12[5], 0)
  expect_equal(c(torus_dist(x = rbind(fit11[1:2]), y = rbind(mu11[1:2])),
                 fit11[3:5] - xi11), rep(0, 4), tolerance = 0.5)
  expect_equal(c(torus_dist(x = rbind(fit12[1:2]), y = rbind(mu12[1:2])),
                 fit12[3:5] - xi12), rep(0, 4), tolerance = 0.5)

})

test_that("Recover parameters under homogeneity and independence via MLE", {

  fit11 <- unname(unlist(fit_bwc_mle(samp11, hom = TRUE, indep = TRUE)[1:5]))
  fit12 <- unname(unlist(fit_bwc_mle(samp12, hom = TRUE, indep = TRUE)[1:5]))

  expect_equal(fit11[5], 0)
  expect_equal(fit12[5], 0)
  expect_equal(c(torus_dist(x = rbind(fit11[1:2]), y = rbind(mu11[1:2])),
                 fit11[3:5] - xi11), rep(0, 4), tolerance = 0.5)
  expect_equal(c(torus_dist(x = rbind(fit12[1:2]), y = rbind(mu12[1:2])),
                 fit12[3:5] - xi12), rep(0, 4), tolerance = 0.5)

})

test_that("BWC fit errors", {

  expect_error(fit_bwc_mle(x = c(0, 0), start = 4))

})
