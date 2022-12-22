
#' @title Fourier expansion of a given curve
#'
#' @description Computation of the Fourier expansion coefficients of a
#' given curve.
#'
#' @param curve points of the curve.
#' @param norm_prop percentage of explained norm. Defaults to \code{1}.
#' @param N number of Gaussian quadrature points, passed to
#' \link[=sphunif]{Gauss_Legen_nodes}. Defaults to \code{1280}.
#' @param K number of terms in the Fourier expansion. Defaults to \code{15}.
#' @param at2 do the \code{atan2} fit instead of the sine fit (only using
#' \eqn{S_m})? Defaults to \code{TRUE}. \code{at2 = FALSE} is not
#' recommended to use.
#' @return The coefficients of the fit (see \code{\link{ridge_curve}}). A list
#' with entries:
#' \item{cos_a}{contains \eqn{a_0,a_1,\ldots,a_m}.}
#' \item{sin_b}{contains \eqn{b_1,\ldots,b_m}.}
#' @examples
#' \donttest{
#' # Zero mean
#' ridge0 <- ridge_bvm(mu = c(0, 0), kappa = c(1, 2, -5), subint_1 = 5e2,
#'                     subint_2 = 5e2)
#' coefs <- ridge_fourier_fit(ridge0)
#' th <- seq(-pi, pi, l = 500)
#' plot(ridge0, xlim = c(-pi, pi), ylim = c(-pi, pi))
#' points(ridge_curve(th, mu = c(0, 0), coefs = coefs, at2 = TRUE), col = 3,
#'        cex = 0.5)
#'
#' # Non-zero mean from a zero-mean ridge
#' mu <- c(1.4, 2)
#' ridge1 <- ridge_bvm(mu = mu, kappa = c(1, 2, -5), subint_1 = 5e2,
#'                     subint_2 = 5e2) # Just for plot
#' plot(ridge1, xlim = c(-pi, pi), ylim = c(-pi, pi))
#' points(mu[1], mu[2], col = 4, pch = "*", cex = 5)
#' points(ridge_curve(th, mu = mu, coefs = coefs), col = 3, cex = 0.5)
#'
#' # Other zero-mean example
#' mu <- c(0, 0)
#' ridge <- ridge_bwc(mu = mu, xi = c(0.3, 0.5, 0.7), subint_1 = 5e2,
#'                    subint_2 = 5e2)
#' plot(ridge, xlim = c(-pi, pi), ylim = c(-pi, pi))
#' coefs <- ridge_fourier_fit(ridge)
#' points(ridge_curve(th, mu = mu, coefs = coefs), col = 4, cex = 0.5)
#'
#' # Another zero-mean example
#' mu <- c(0, 0)
#' ridge <- ridge_bwc(mu = mu, xi = c(0.8, 0.1, 0.75), subint_1 = 5e2,
#'                    subint_2 = 5e2)
#' plot(ridge, xlim = c(-pi, pi), ylim = c(-pi, pi))
#' coefs <- ridge_fourier_fit(ridge)
#' points(ridge_curve(th, mu = mu, coefs = coefs), col = 4, cex = 0.5)}
#' @export
ridge_fourier_fit <- function(curve, K = 15, norm_prop = 1, N = 1280,
                              at2 = TRUE) {

  # Calculate nodes and weights
  nodes <- drop(sphunif::Gauss_Legen_nodes(a = -pi, b = pi, N = N))
  w <- drop(sphunif::Gauss_Legen_weights(a = -pi, b = pi, N = N))

  # Find the closest point of the curve to each node and save its y value
  y <- rep(0, length = N)
  for (i in 1:N) {

    y[i] <- curve[which.min(abs(curve[, 1] - nodes[i])), 2]

  }

  # Fourier atan2 or sine fit?
  if (at2) {

    # Estimate cos/sin coefficients
    y_cos_w <- cos(y) * w / pi
    y_sin_w <- sin(y) * w / pi
    cos_a <- sapply(0:K, function(k) sum(cos(k * nodes) * y_cos_w))
    sin_b <- sapply(1:K, function(k) sum(sin(k * nodes) * y_sin_w))

    # Squared norm of the curve and norm cutoff
    norm_sq_cos <- pi * (cos_a[1]^2 / 2 + cumsum(cos_a[-1]^2))
    norm_sq_sin <- pi * cumsum(sin_b^2)
    cutoff_cos <- which(norm_sq_cos / norm_sq_cos[K] >= norm_prop)
    cutoff_sin <- which(norm_sq_sin / norm_sq_sin[K] >= norm_prop)
    cutoff <- max(cutoff_cos, cutoff_sin)

    # Coefficients
    coefs <- list("cos_a" = cos_a[1:(cutoff + 1)], "sin_b" = sin_b[1:cutoff])

  } else {

    # Estimate sine coefficients
    y_sin_w <- sin(y) * w / pi
    coefs <- sapply(1:K, function(k) sum(sin(k * nodes) * y_sin_w))

    # Squared norm of the curve and norm cutoff
    norm_sq <- pi * cumsum(coefs^2)
    cutoff <- which(norm_sq / norm_sq[K] >= norm_prop)
    coefs <- coefs[1:cutoff]

  }

  # Coefficients
  return(coefs)

}
