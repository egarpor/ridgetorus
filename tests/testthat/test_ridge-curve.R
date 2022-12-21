
test_that("der_ridge_curve does the correct derivative for at2 = FALSE", {

  for (i in 1:50) {

    set.seed(i)
    th <- runif(1, -pi, pi)
    mu <- runif(2, -pi, pi)
    coefs <- rnorm(rpois(1, lambda = 3) + 1)
    ind_var <- sample(1:2, size = 1)
    expect_equal(
      drop(numDeriv::jacobian(func = ridge_curve, x = th, mu = mu,
                              coefs = coefs, ind_var = ind_var, at2 = FALSE)),
      drop(der_ridge_curve(theta = th, mu = mu, coefs = coefs,
                           ind_var = ind_var, at2 = FALSE)),
      tolerance = 1e-4
    )

  }

})

test_that("dist_ridge_curve equals for der = TRUE/FALSE for at2 = FALSE ", {

  for (i in 1:50) {

    set.seed(i)
    a <- runif(2, -pi, pi)
    mu <- runif(2, -pi, pi)
    coefs <- rnorm(rpois(1, lambda = 3) + 1)
    ind_var <- sample(1:2, size = 1)
    expect_equal(
      dist_ridge_curve(a, mu = mu, coefs = coefs, ind_var = ind_var,
                       N = 5e2, der = TRUE, at2 = FALSE),
      dist_ridge_curve(a, mu = mu, coefs = coefs, ind_var = ind_var,
                       N = 5e2, der = FALSE, at2 = FALSE),
      tolerance = 1e-3
    )

  }

})

test_that("arclength_ridge_curve rectifies curves for at2 = FALSE", {

  for (i in 1:20) {

    set.seed(i)
    mu <- runif(2, -pi, pi)
    coefs <- rnorm(rpois(1, lambda = 2) + 1)
    ind_var <- sample(1:2, size = 1)
    n <- 100
    alpha <- arclength_ridge_curve(mu = mu, coefs = coefs, ind_var = ind_var,
                                   N = n, at2 = FALSE)
    expect_equal(drop(ridge_curve(theta = alpha[1], mu = mu, coefs = coefs,
                                  ind_var = ind_var, at2 = FALSE)), mu)
    ind <- sample(n, size = 2)
    max_length <- dist_ridge_curve(alpha = c(0, 2 * pi), mu = mu, coefs = coefs,
                                   ind_var = ind_var, N = 5e2,
                                   shortest = FALSE, at2 = FALSE)
    expect_equal(dist_ridge_curve(alpha = alpha[ind], mu = mu, coefs = coefs,
                                  ind_var = ind_var, N = 5e2,
                                  shortest = FALSE, at2 = FALSE) / max_length,
                 abs(diff(ind)) / n, tolerance = 1e-2)

  }

})

test_that("der_ridge_curve does the correct derivative", {

  for (i in 1:50) {

    set.seed(i)
    th <- runif(1, -pi, pi)
    mu <- runif(2, -pi, pi)
    K <- rpois(1, lambda = 3) + 1
    coefs <- list("cos_a" = rnorm(K + 1), "sin_b" = rnorm(K))
    ind_var <- sample(1:2, size = 1)
    expect_equal(
      drop(numDeriv::jacobian(func = ridge_curve, x = th, mu = mu,
                              coefs = coefs, ind_var = ind_var)),
      drop(der_ridge_curve(theta = th, mu = mu, coefs = coefs,
                           ind_var = ind_var))
    )

  }

})

test_that("der_ridge_curve normalizes correctly", {

  for (i in 1:50) {

    set.seed(i)
    th <- runif(1, -pi, pi)
    mu <- runif(2, -pi, pi)
    K <- rpois(1, lambda = 3) + 1
    coefs <- list("cos_a" = rnorm(K + 1), "sin_b" = rnorm(K))
    ind_var <- sample(1:2, size = 1)
    der <- der_ridge_curve(theta = th, mu = mu, coefs = coefs,
                           ind_var = ind_var, norm = 1)
    expect_equal(sqrt(rowSums(der^2)), 1)

  }

})

test_that("dist_ridge_curve equals for der = TRUE/FALSE", {

  for (i in 1:2) {

    set.seed(i)
    a <- runif(2, -pi, pi)
    mu <- runif(2, -pi, pi)
    K <- rpois(1, lambda = 3) + 1
    coefs <- list("cos_a" = rnorm(K + 1), "sin_b" = rnorm(K))
    ind_var <- sample(1:2, size = 1)
    expect_equal(
      dist_ridge_curve(a, mu = mu, coefs = coefs, ind_var = ind_var,
                       N = 5e2, der = TRUE),
      dist_ridge_curve(a, mu = mu, coefs = coefs, ind_var = ind_var,
                       N = 5e2, der = FALSE),
      tolerance = 1e-3
    )

  }

})

test_that("arclength_ridge_curve rectifies curves", {

  for (i in 1:20) {

    set.seed(i)
    mu <- runif(2, -pi, pi)
    K <- rpois(1, lambda = 3) + 1
    coefs <- list("cos_a" = rnorm(K + 1), "sin_b" = rnorm(K))
    ind_var <- sample(1:2, size = 1)
    n <- 100
    alpha <- arclength_ridge_curve(mu = mu, coefs = coefs, ind_var = ind_var,
                                    N = n)
    expect_equal(drop(ridge_curve(theta = alpha[1], mu = mu, coefs = coefs,
                                  ind_var = ind_var)), mu)
    ind <- sample(n, size = 2)
    max_length <- dist_ridge_curve(alpha = c(0, 2 * pi), mu = mu, coefs = coefs,
                                   ind_var = ind_var, N = 5e2,
                                   shortest = FALSE)
    expect_equal(dist_ridge_curve(alpha = alpha[ind], mu = mu, coefs = coefs,
                                  ind_var = ind_var, N = 5e2,
                                  shortest = FALSE) / max_length,
                 abs(diff(ind)) / n, tolerance = 2e-1)

  }

})
