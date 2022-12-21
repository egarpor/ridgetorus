
# Check symmetry properties of the ridge

## Translation

# Bivariate von Mises
kcase1 <- c(0.8, 0.45, 0.9)
mu1 <- c(1, 2)
ridge1 <- ridge_bvm(mu = c(0, 0), kappa = kcase1, subint_1 = 5e2,
                    subint_2 = 5e2)
ridge1b <- ridge_bvm(mu = mu1, kappa = kcase1, subint_1 = 5e2, subint_2 = 5e2)
ridge1b[, 1] <- sdetorus::toPiInt(ridge1b[, 1] - mu1[1])
ridge1b[, 2] <- sdetorus::toPiInt(ridge1b[, 2] - mu1[2])

kcase2 <- c(0.1, 0.2, 0.4)
mu2 <- c(-2, 0.5)
ridge2 <- ridge_bvm(mu = c(0, 0), kappa = kcase2, subint_1 = 5e2,
                    subint_2 = 5e2)
ridge2b <- ridge_bvm(mu = mu2, kappa = kcase2, subint_1 = 5e2, subint_2 = 5e2)
ridge2b[, 1] <- sdetorus::toPiInt(ridge2b[, 1] - mu2[1])
ridge2b[, 2] <- sdetorus::toPiInt(ridge2b[, 2] - mu2[2])

kcase3 <- c(1.5, 1.1, 1.2)
mu3 <- c(3, -1)
ridge3 <- ridge_bvm(mu = c(0, 0), kappa = kcase3, subint_1 = 5e2,
                    subint_2 = 5e2)
ridge3b <- ridge_bvm(mu = mu3, kappa = kcase3, subint_1 = 5e2, subint_2 = 5e2)
ridge3b[, 1] <- sdetorus::toPiInt(ridge3b[, 1] - mu3[1])
ridge3b[, 2] <- sdetorus::toPiInt(ridge3b[, 2] - mu3[2])

test_that("Ridge is invariant under translation", {

  expect_equal(max(torus_dist(ridge1, ridge1b)), 0, tolerance = 0.01)
  expect_equal(max(torus_dist(ridge2, ridge2b)), 0, tolerance = 0.01)
  expect_equal(max(torus_dist(ridge3, ridge3b)), 0, tolerance = 0.01)

})

# Bivariate wrapped Cauchy
xicase1 <- c(0.5, 0.3, 0.7)
mu1 <- c(0.5, 1.4)
ridge1 <- ridge_bwc(mu = c(0, 0), xi = xicase1, subint_1 = 5e2, subint_2 = 5e2)
ridge1b <- ridge_bwc(mu = mu1, xi = xicase1, subint_1 = 5e2, subint_2 = 5e2)
ridge1b[, 1] <- sdetorus::toPiInt(ridge1b[, 1] - mu1[1])
ridge1b[, 2] <- sdetorus::toPiInt(ridge1b[, 2] - mu1[2])

xicase2 <- c(0.1, 0.2, 0.4)
mu2 <- c(-2.3, 0.95)
ridge2 <- ridge_bwc(mu = c(0, 0), xi = xicase2, subint_1 = 5e2, subint_2 = 5e2)
ridge2b <- ridge_bwc(mu = mu2, xi = xicase2, subint_1 = 5e2, subint_2 = 5e2)
ridge2b[, 1] <- sdetorus::toPiInt(ridge2b[, 1] - mu2[1])
ridge2b[, 2] <- sdetorus::toPiInt(ridge2b[, 2] - mu2[2])

xicase3 <- c(0.7, 0.6, 0.25)
mu3 <- c(3, -0.1)
ridge3 <- ridge_bwc(mu = c(0, 0), xi = xicase3, subint_1 = 5e2, subint_2 = 5e2)
ridge3b <- ridge_bwc(mu = mu3, xi = xicase3, subint_1 = 5e2, subint_2 = 5e2)
ridge3b[, 1] <- sdetorus::toPiInt(ridge3b[, 1] - mu3[1])
ridge3b[, 2] <- sdetorus::toPiInt(ridge3b[, 2] - mu3[2])

test_that("Ridge is invariant under translation", {

  expect_equal(max(torus_dist(ridge1, ridge1b)), 0, tolerance = 0.01)
  expect_equal(max(torus_dist(ridge2, ridge2b)), 0, tolerance = 0.01)
  expect_equal(max(torus_dist(ridge3, ridge3b)), 0, tolerance = 0.01)

})

## Rotation

# Bivariate von Mises
kcase1 <- c(0.8, 0.45, 0.9)
kcase12 <- c(0.45, 0.8, 0.9)
ridge1 <- ridge_bvm(mu = c(0, 0), kappa = kcase1, subint_1 = 5e2,
                    subint_2 = 5e2)
ridge1b <- ridge_bvm(mu = c(0, 0), kappa = kcase12, subint_1 = 5e2,
                     subint_2 = 5e2)

kcase2 <- c(0.1, 0.2, 0.4)
kcase22 <- c(0.2, 0.1, 0.4)
ridge2 <- ridge_bvm(mu = c(0, 0), kappa = kcase2, subint_1 = 5e2,
                    subint_2 = 5e2)
ridge2b <- ridge_bvm(mu = c(0, 0), kappa = kcase22, subint_1 = 5e2,
                     subint_2 = 5e2)

kcase3 <- c(1.5, 1.1, 1.2)
kcase32 <- c(1.1, 1.5, 1.2)
ridge3 <- ridge_bvm(mu = c(0, 0), kappa = kcase3, subint_1 = 5e2,
                    subint_2 = 5e2)
ridge3b <- ridge_bvm(mu = c(0, 0), kappa = kcase32, subint_1 = 5e2,
                     subint_2 = 5e2)

test_that("Ridge is symmetric under rotation", {

  expect_equal(max(torus_dist(ridge1, ridge1b[, 2:1])), 0, tolerance = 0.01)
  expect_equal(max(torus_dist(ridge2, ridge2b[, 2:1])), 0, tolerance = 0.01)
  expect_equal(max(torus_dist(ridge3, ridge3b[, 2:1])), 0, tolerance = 0.01)

})

# Bivariate wrapped Cauchy
xicase1 <- c(0.5, 0.3, 0.7)
xicase12 <- c(0.3, 0.5, 0.7)
ridge1 <- ridge_bwc(mu = c(0, 0), xi = xicase1, subint_1 = 5e2, subint_2 = 5e2)
ridge1b <- ridge_bwc(mu = c(0, 0), xi = xicase12, subint_1 = 5e2,
                     subint_2 = 5e2)

xicase2 <- c(0.1, 0.2, 0.4)
xicase22 <- c(0.2, 0.1, 0.4)
ridge2 <- ridge_bwc(mu = c(0, 0), xi = xicase2, subint_1 = 5e2, subint_2 = 5e2)
ridge2b <- ridge_bwc(mu = c(0, 0), xi = xicase22, subint_1 = 5e2,
                     subint_2 = 5e2)

xicase3 <- c(0.7, 0.6, 0.25)
xicase32 <- c(0.6, 0.7, 0.25)
ridge3 <- ridge_bwc(mu = c(0, 0), xi = xicase3, subint_1 = 5e2, subint_2 = 5e2)
ridge3b <- ridge_bwc(mu = c(0, 0), xi = xicase32, subint_1 = 5e2,
                     subint_2 = 5e2)

test_that("Ridge is symmetric under rotation", {

  expect_equal(max(torus_dist(ridge1, ridge1b[, 2:1])), 0, tolerance = 0.01)
  expect_equal(max(torus_dist(ridge2, ridge2b[, 2:1])), 0, tolerance = 0.01)
  expect_equal(max(torus_dist(ridge3, ridge3b[, 2:1])), 0, tolerance = 0.01)

})

## Reflection

# Bivariate von Mises
kcase1 <- c(0.8, 0.45, 0.8)
kcase12 <- c(0.8, 0.45, -0.8)
ridge1 <- ridge_bvm(mu = c(0, 0), kappa = kcase1, subint_1 = 5e2,
                    subint_2 = 5e2)
ridge1b <- ridge_bvm(mu = c(0, 0), kappa = kcase12, subint_1 = 5e2,
                     subint_2 = 5e2)

kcase2 <- c(0.1, 0.2, 0.4)
kcase22 <- c(0.1, 0.2, -0.4)
ridge2 <- ridge_bvm(mu = c(0, 0), kappa = kcase2, subint_1 = 5e2,
                    subint_2 = 5e2)
ridge2b <- ridge_bvm(mu = c(0, 0), kappa = kcase22, subint_1 = 5e2,
                     subint_2 = 5e2)

kcase3 <- c(1.5, 1.1, 1.2)
kcase32 <- c(1.5, 1.1, -1.2)
ridge3 <- ridge_bvm(mu = c(0, 0), kappa = kcase3, subint_1 = 5e2,
                    subint_2 = 5e2)
ridge3b <- ridge_bvm(mu = c(0, 0), kappa = kcase32, subint_1 = 5e2,
                     subint_2 = 5e2)

test_that("Ridge is symmetric under reflection", {

  # Negative sign on the lowest kappa variable
  expect_equal(max(torus_dist(ridge1, cbind(ridge1b[, 1], -ridge1b[, 2]))), 0,
               tolerance = 0.01)
  expect_equal(max(torus_dist(ridge2, cbind(-ridge2b[, 1], ridge2b[, 2]))), 0,
               tolerance = 0.01)
  expect_equal(max(torus_dist(ridge3, cbind(ridge3b[, 1], -ridge3b[, 2]))), 0,
               tolerance = 0.01)

})

# Bivariate wrapped Cauchy
xicase1 <- c(0.5, 0.3, 0.7)
xicase12 <- c(0.5, 0.3, -0.7)
ridge1 <- ridge_bwc(mu = c(0, 0), xi = xicase1, subint_1 = 5e2, subint_2 = 5e2)
ridge1b <- ridge_bwc(mu = c(0, 0), xi = xicase12, subint_1 = 5e2,
                     subint_2 = 5e2)

xicase2 <- c(0.1, 0.25, 0.4)
xicase22 <- c(0.1, 0.25, -0.4)
ridge2 <- ridge_bwc(mu = c(0, 0), xi = xicase2, subint_1 = 5e2, subint_2 = 5e2)
ridge2b <- ridge_bwc(mu = c(0, 0), xi = xicase22, subint_1 = 5e2,
                     subint_2 = 5e2)

xicase3 <- c(0.75, 0.6, 0.25)
xicase32 <- c(0.75, 0.6, -0.25)
ridge3 <- ridge_bwc(mu = c(0, 0), xi = xicase3, subint_1 = 5e2, subint_2 = 5e2)
ridge3b <- ridge_bwc(mu = c(0, 0), xi = xicase32, subint_1 = 5e2,
                     subint_2 = 5e2)

test_that("Ridge is symmetric under reflection", {

  # Negative sign on the lowest kappa variable
  expect_equal(max(torus_dist(ridge1, cbind(ridge1b[, 1], -ridge1b[, 2]))), 0,
               tolerance = 0.01)
  expect_equal(max(torus_dist(ridge2, cbind(-ridge2b[, 1], ridge2b[, 2]))), 0,
               tolerance = 0.01)
  expect_equal(max(torus_dist(ridge3, cbind(ridge3b[, 1], -ridge3b[, 2]))), 0,
               tolerance = 0.01)

})

## Check that the implicit equation is 0

# Bivariate von Mises
mu1 <- c(-pi, 1)
kappa1 <- c(0.3, 1, -0.2)
mu2 <- c(1, 2)
kappa2 <- c(0.4, 0.1, 0.8)
kappa32 <- c(1, 1, 0.5)

# Normal cases
ridge1 <- ridge_bvm(mu = mu1, kappa = kappa1, subint_1 = 5e2, subint_2 = 5e2)
ridge2 <- ridge_bvm(mu = mu2, kappa = kappa2, subint_1 = 5e2, subint_2 = 5e2)

# With evaluation points
ridge3 <- ridge_bvm(mu = mu2, kappa = kappa2, subint_1 = 5e2, subint_2 = 5e2,
                    eval_points = seq(-2, 2, l = 100))
ridge32 <- ridge_bvm(mu = mu2, kappa = kappa32, subint_1 = 5e2, subint_2 = 5e2,
                    eval_points = seq(-2, 2, l = 100))

# With limit cases
mu4 <- c(0, 0)
kappa4 <- c(0.5, 0.5, 0.8)
ridge4 <- ridge_bvm(mu = mu4, kappa = kappa4, subint_1 = 5e2, subint_2 = 5e2)
mu5 <- c(0, 0)
kappa5 <- c(0.5, 0.5, 0.3)
ridge5 <- ridge_bvm(mu = mu5, kappa = kappa5, subint_1 = 5e2, subint_2 = 5e2)

val1 <- apply(sdetorus::toPiInt(t(t(ridge1) - mu1)), 1, function(x)
  ridgetorus:::implicit_equation(
    theta2 = x[2], theta1 = x[1], kappa = kappa1[1:2],
    Lambda = matrix(c(0, kappa1[3], kappa1[3], 0), nrow = 2), density = "bvm"))
val2 <- apply(sdetorus::toPiInt(t(t(ridge2) - mu2)), 1, function(x)
  ridgetorus:::implicit_equation(
    theta2 = x[2], theta1 = x[1], kappa = kappa2[1:2],
    Lambda = matrix(c(0, kappa2[3], kappa2[3], 0), nrow = 2), density = "bvm"))
val3 <- apply(sdetorus::toPiInt(t(t(ridge3) - mu2)), 1, function(x)
  ridgetorus:::implicit_equation(
    theta2 = x[2], theta1 = x[1], kappa = kappa2[1:2],
    Lambda = matrix(c(0, kappa2[3], kappa2[3], 0), nrow = 2), density = "bvm"))
val32 <- apply(sdetorus::toPiInt(t(t(ridge32) - mu2)), 1, function(x)
  ridgetorus:::implicit_equation(
    theta2 = x[2], theta1 = x[1], kappa = kappa32[1:2],
    Lambda = matrix(c(0, kappa32[3], kappa32[3], 0), nrow = 2),
    density = "bvm"))
val4 <- apply(sdetorus::toPiInt(t(t(ridge4) - mu4)), 1, function(x)
  ridgetorus:::implicit_equation(
    theta2 = x[2], theta1 = x[1], kappa = kappa4[1:2],
    Lambda = matrix(c(0, kappa4[3], kappa4[3], 0), nrow = 2), density = "bvm"))
val5 <- apply(sdetorus::toPiInt(t(t(ridge5) - mu5)), 1, function(x)
  ridgetorus:::implicit_equation(
    theta2 = x[2], theta1 = x[1], kappa = kappa5[1:2],
    Lambda = matrix(c(0, kappa5[3], kappa5[3], 0), nrow = 2), density = "bvm"))

test_that("Implicit is 0", {

  expect_equal(val1, rep(0, length(ridge1[, 1])), tolerance = 0.01)
  expect_equal(val2, rep(0, length(ridge2[, 1])), tolerance = 0.01)
  expect_equal(val3, rep(0, length(ridge3[, 1])), tolerance = 0.01)
  expect_equal(val32, rep(0, length(ridge32[, 1])), tolerance = 0.01)
  expect_equal(val4, rep(0, length(ridge4[, 1])), tolerance = 0.01)
  expect_equal(val5, rep(0, length(ridge5[, 1])), tolerance = 0.01)

})

# Bivariate wrapped Cauchy
mu1 <- c(-pi, 1)
xi1 <- c(0.3, 0.5, 0.2)
mu2 <- c(1, 2)
xi2 <- c(0.4, 0.1, -0.8)
mu3 <- c(0, 0)
xi3 <- c(0.3, 0.3, 0.5)
ridge1 <- ridge_bwc(mu = mu1, xi = xi1,  subint_1 = 5e2, subint_2 = 5e2)
ridge2 <- ridge_bwc(mu = mu2, xi = xi2, subint_1 = 5e2, subint_2 = 5e2)
ridge3 <- ridge_bwc(mu = mu2, xi = xi2, subint_1 = 5e2, subint_2 = 5e2,
                    eval_points = seq.int(-3, 3))
ridge4 <- ridge_bwc(mu = mu3, xi = xi3, subint_1 = 5e2, subint_2 = 5e2)
val1 <- apply(sdetorus::toPiInt(t(t(ridge1) - mu1)), 1, function(x)
  ridgetorus:::implicit_equation(theta2 = x[2], theta1 = x[1], xi = xi1,
                                 density = "bwc"))
val2 <- apply(sdetorus::toPiInt(t(t(ridge2) - mu2)), 1, function(x)
  ridgetorus:::implicit_equation(theta2 = x[2], theta1 = x[1], xi = xi2,
                                 density = "bwc"))
val3 <- apply(sdetorus::toPiInt(t(t(ridge3) - mu2)), 1, function(x)
  ridgetorus:::implicit_equation(theta2 = x[2], theta1 = x[1], xi = xi2,
                                 density = "bwc"))
val4 <- apply(sdetorus::toPiInt(t(t(ridge4) - mu3)), 1, function(x)
  ridgetorus:::implicit_equation(theta2 = x[2], theta1 = x[1], xi = xi3,
                                 density = "bwc"))

test_that("Implicit is 0", {

  expect_equal(val1, rep(0, length(ridge1[, 1])), tolerance = 0.01)
  expect_equal(val2, rep(0, length(ridge2[, 1])), tolerance = 0.01)
  expect_equal(val3, rep(0, length(ridge3[, 1])), tolerance = 0.01)
  expect_equal(val4, rep(0, length(ridge4[, 1])), tolerance = 0.01)

})

# Bivariate wrapped normal
mu1 <- c(0, 0)
Sigma1 <- matrix(data = c(10, 2, 2, 5), nrow = 2)
mu2 <- c(0, 0)
Sigma2 <- matrix(data = c(2, 1, 1, 1), nrow = 2)
k <- c(-1, 0, 1)
K <- as.matrix(expand.grid(k, k))

ridge1 <- ridge_bwn(mu = mu1, Sigma = Sigma1, subint_1 = 5e2, subint_2 = 5e2)
ridge2 <- ridge_bwn(mu = mu2, Sigma = Sigma2, subint_1 = 5e2, subint_2 = 5e2)
ridge3 <- ridge_bwn(mu = mu2, Sigma = Sigma2, subint_1 = 5e2, subint_2 = 5e2,
                    eval_points = seq.int(-3, 3))

val1 <- apply(sdetorus::toPiInt(t(t(ridge1))), 1,
              function(x)  ridgetorus:::implicit_equation(
                theta2 = x[2], theta1 = x[1], Sigma = Sigma1, k = K, mu = mu1,
                density = "bwn"))
val2 <- apply(sdetorus::toPiInt(t(t(ridge2))), 1,
              function(x) ridgetorus:::implicit_equation(
                theta2 = x[2], theta1 = x[1], Sigma = Sigma2, k = K, mu = mu2,
                density = "bwn"))
val3 <- apply(sdetorus::toPiInt(t(t(ridge3))), 1,
              function(x) ridgetorus:::implicit_equation(
                theta2 = x[2], theta1 = x[1], Sigma = Sigma2, k = K, mu = mu2,
                density = "bwn"))

test_that("Implicit is 0", {

  expect_equal(val1, rep(0, length(ridge1[, 1])), tolerance = 0.01)
  expect_equal(val2, rep(0, length(ridge2[, 1])), tolerance = 0.01)
  expect_equal(val3, rep(0, length(ridge3[, 1])), tolerance = 0.01)

})

# Errors and warnings tests
test_that("BWC errors", {

  expect_error(ridge_bwc(mu = c(0, 0), xi = c(0.3, 0.4, 0.5), subint_1 = -1,
                         subint_2 = 10))
  expect_error(ridge_bwc(mu = c(0, 0), xi = c(0.3, 0.4, 0.5), subint_1 = 10,
                         subint_2 = -10))
  expect_error(ridge_bwc(mu = c(0, 0), xi = c(0.3, 0.4, 2.5), subint_1 = 5e2,
                         subint_2 = 5e2))
  expect_error(ridge_bwc(mu = c(0, 0), xi = c(0.3, 1.4, 0.5), subint_1 = 5e2,
                         subint_2 = 5e2))
  expect_error(ridge_bwc(mu = c(0, 0), xi = c(1.3, 0.4, 0.5), subint_1 = 5e2,
                         subint_2 = 5e2))
  expect_error(ridge_bwc(mu = c(0, 0), xi = c(0.3, 0.4, 0.5, 0), subint_1 = 5e2,
                         subint_2 = 5e2))
  expect_error(ridge_bwc(mu = c(1), xi = c(0.3, 0.4, 0.5), subint_1 = 5e2,
                         subint_2 = 5e2))

})

test_that("BvM errors", {

  expect_error(ridge_bvm(mu = c(0, 0), kappa = c(0.3, 0.4, 0), subint_1 = -1,
                         subint_2 = 10))
  expect_error(ridge_bvm(mu = c(0, 0), kappa = c(0.3, 0.4, 0), subint_1 = 1,
                         subint_2 = -10))
  expect_error(ridge_bvm(mu = c(0, 0), kappa = c(-0.3, 0.4, 2.5), subint_1 = 1,
                         subint_2 = 1))
  expect_error(ridge_bvm(mu = c(0, 0), kappa = c(0.3, -1.4, 0.5), subint_1 = 1,
                         subint_2 = 1))
  expect_error(ridge_bvm(mu = c(0, 0), kappa = c(1.3, 0, 0, 0), subint_1 = 5e2,
                         subint_2 = 5e2))
  expect_error(ridge_bvm(mu = c(1), kappa = c(0.3, 0.4, 0.5), subint_1 = 5e2,
                         subint_2 = 5e2))

})

test_that("BWN errors", {

  expect_error(ridge_bwn(mu = c(0, 0), subint_1 = -1, subint_2 = 10,
                         Sigma = c(10, 2, 2, 5)))

  expect_error(ridge_bwn(mu = c(0, 0), subint_1 = -1, subint_2 = 10,
                         Sigma = matrix(data = c(10, 2, 2, 5, 0, 0), nrow = 3)
  ))
  expect_error(ridge_bwn(mu = c(0, 0, 0), subint_1 = -1, subint_2 = 10,
                         Sigma = matrix(data = c(10, 2, 2, 5), nrow = 2)
  ))
  expect_error(ridge_bwn(mu = c(0, 0), subint_1 = -1, subint_2 = 10,
                         Sigma = matrix(data = c(10, 2, 2, 5), nrow = 2)))

  expect_error(ridge_bwn(mu = c(0, 0), subint_1 = 100, subint_2 = -10,
                         Sigma = matrix(data = c(10, 2, 2, 5), nrow = 2)))

})

test_that("implicit_equation() errors", {

  expect_error(implicit_equation(theta2 = 0, theta1 = 0, density = "model"))

})
