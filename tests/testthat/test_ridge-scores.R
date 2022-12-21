
test_that("torus_dist() fails with incorrect arguments", {

  expect_error(torus_dist(x = 1:2, y = rbind(4:5, 4:5)))
  expect_error(torus_dist(x = rbind(1:3, 3:1), y = 4:5))
  expect_error(torus_dist(x = rbind(1:3, 3:1), y = rbind(4:5, 4:5)))

})

test_that("torus_dist() basic checks", {

  expect_equal(torus_dist(x = rbind(1:2, 2:1), y = 4:5),
               c(torus_dist(x = rbind(1:2), y = 4:5),
                 torus_dist(x = rbind(2:1), y = 4:5)))
  expect_equal(torus_dist(x = rbind(1:2, 3:2), y = rbind(4:5, 5:4)),
               c(torus_dist(x = rbind(1:2), y = 4:5),
                 torus_dist(x = rbind(3:2), y = 5:4)))

})

test_that("first scores are scaled", {

  mu <- c(-0.5, 1.65)
  th <- seq(-pi, pi, l = 200)
  K <- 5
  coefs <- list(cos_a = 1 / (1:(K + 1))^3, sin_b = 1 / (1:K)^3)
  n <- 10
  x <- ridge_curve(theta = th, mu = mu, coefs = coefs)
  s_uns <- max(ridge_scores(x, mu = mu, coefs = coefs, scale = FALSE)$scores)
  s_sc <- max(ridge_scores(x, mu = mu, coefs = coefs, scale = TRUE)$scores)
  expect_equal(max(max(s_uns), pi), s_uns)
  expect_equal(max(s_sc, pi), pi)

})

test_that("max_score_2 is invariant to mu", {

  mu <- runif(2, min = -pi, pi)
  coefs <- list(cos_a = runif(2), sin_b = runif(1))
  for (ind_var in 1:2) {
    expect_equal(max_score_2(mu = mu, coefs = coefs, ind_var = ind_var),
                 max_score_2(mu = mu + runif(2), coefs = coefs,
                             ind_var = ind_var),
                 tolerance = 1e-2)
  }

})

test_that("max_score_2 gives maximal score", {

  mu <- runif(2, min = -pi, pi)
  coefs <- list(cos_a = runif(2), sin_b = runif(1))
  th <- seq(-pi, pi, l = 101)[-101]
  th <- as.matrix(expand.grid(th, th))
  projs <- proj_ridge_curve(x = th, mu = mu, coefs = coefs, ind_var = 1)
  projs <- ridge_curve(theta = projs$theta_proj, mu = mu, coefs = coefs,
                       ind_var = 1)
  expect_equal(max(torus_dist(x = th, y = projs)),
               max_score_2(mu = mu, coefs = coefs, ind_var = 1, L = 25, f = 2),
               tolerance = 5e-2)
  i <- which.max(torus_dist(x = th, y = projs))
  projs2 <- proj_ridge_curve(x = th, mu = mu, coefs = coefs, ind_var = 1,
                            arclength = TRUE)
  projs2 <- ridge_curve(theta = projs2$theta_proj, mu = mu, coefs = coefs,
                       ind_var = 1)
  expect_equal(max(torus_dist(x = th, y = projs2)),
               max_score_2(mu = mu, coefs = coefs, ind_var = 1, L = 25, f = 2),
               tolerance = 5e-2)
  i <- which.max(torus_dist(x = th, y = projs2))

})
