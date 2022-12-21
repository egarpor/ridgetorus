
# Check that atan2 fit reproduces the ridge

K <- 15
eps <- 1e-2

## Bivariate von Mises

# Second concentration has to be larger as ind_var is hardcoded
kappa <- rbind(c(1, 2, 0),
               c(2, 3, 5),
               c(1, 2, 10),
               c(0.2, 0.3, 0.5),
               c(0.3, 0.95, -0.5),
               c(0.1, 0.7, -0.8),
               c(0.01, 0.4, 0.02))

for (i in seq_len(nrow(kappa))) {

  ridge <- ridge_bvm(mu = c(0, 0), kappa = kappa[i, ],
                     subint_1 = 5e2, subint_2 = 5e2)
  fit <- ridge_fourier_fit(curve = ridge, at2 = TRUE, K = K)
  coefs <- list(cos_a = fit$cos_a, sin_b = fit$sin_b)
  ridge2 <- ridge_curve(theta = ridge[, 1], coefs = coefs, at2 = TRUE,
                        ind_var = 1)
  max_dist <- max(apply(X = ridge, FUN = function(y)
    min(torus_dist(x = ridge2, y = y)), MARGIN = 1))
  plot(ridge, main = paste("BvM i =", i, "d =", max_dist), cex = 0.5,
       xlim = c(-pi, pi), ylim = c(-pi, pi))
  points(ridge2, col = 4, cex = 0.2, pch = 16)
  test_that(paste("Fit BvM", i, "reproduces the ridge"), {
    expect_true(max_dist < eps)
  })

}

## Bivariate wrapped Cauchy

# Second concentration has to be larger as ind_var is hardcoded
xi <- rbind(c(0.2, 0.4, 0),
            c(0.25, 0.5, 0.25),
            c(0.3, 0.95, -0.5),
            c(0.15, 0.6, -0.9),
            c(0.05, 0.01, 0.9),
            c(0.005, 0.2, 0.01),
            c(0.8, 0.81, 0.85))

for (i in seq_len(nrow(xi))) {

  ridge <- ridge_bwc(mu = c(0, 0), xi = xi[i, ], subint_1 = 5e2, subint_2 = 5e2)
  fit <- ridge_fourier_fit(curve = ridge, at2 = TRUE, K = K)
  coefs <- list(cos_a = fit$cos_a, sin_b = fit$sin_b)
  ridge2 <- ridge_curve(theta = ridge[, 1], coefs = coefs, ind_var = 1,
                        at2 = TRUE)
  max_dist <- max(apply(X = ridge, FUN = function(y)
    min(torus_dist(x = ridge2, y = y)), MARGIN = 1))
  plot(ridge, main = paste("BWC i =", i, "d =", max_dist), cex = 0.5,
       xlim = c(-pi, pi), ylim = c(-pi, pi))
  points(ridge2, col = 4, cex = 0.2, pch = 16)
  test_that(paste("Fit BWC", i, "reproduces the ridge"), {
    expect_true(max_dist < eps)
  })

}
