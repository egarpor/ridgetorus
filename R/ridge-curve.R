
#' @title Fourier-fitted ridge curve and related utilities
#'
#' @description Given the angles \code{theta} in \eqn{[-\pi, \pi)},
#' \code{ridge_curve} computes the Fourier-fitted ridge curve
#' \eqn{(\theta, r_1(\theta))} or \eqn{(r_2(\theta), \theta)}, where
#' \deqn{r_j(\theta):=\mathrm{atan2}(S_m (\theta),
#' C_m (\theta))} with \eqn{C_m(x) :=
#' a_0/2 + \sum_{k=1}^m a_k \cos(kx)} and \eqn{S_m(x) :=
#' \sum_{k=1}^m b_k \sin(kx)} for \eqn{j = 1,2}. \code{der_ridge_curve} and
#' \code{dist_ridge_curve} compute the derivatives of and the distances along
#' these curves, respectively. \code{alpha_ridge_curve} provides a uniform
#' grid of the ridge curve using the arc-length parametrization.
#' \code{proj_ridge_curve} gives the ridge's \eqn{\theta} for which the curve
#' is closer to any point on \eqn{[-\pi, \pi)^2}.
#'
#' @param theta vector \eqn{\theta} of size \code{nth}.
#' @param mu a vector of size \code{2} giving \eqn{(\mu_1, \mu_2)}. Defaults
#' to \code{c(0, 0)}.
#' @param coefs list of coefficients \code{cos_a} (\eqn{a_k}) and
#' \code{sin_b} (\eqn{b_k} giving the Fourier fit of the ridge curve.
#' Defaults to \code{list(cos_a = c(0, 0), sin_b = 0)}. See examples.
#' @param ind_var index \eqn{j} of the variable that parametrizes the ridge.
#' Defaults to \code{1}.
#' @param norm normalize tangent vectors? If different from \code{NULL}
#' (the default), the vectors are normalized to have the given \code{norm}.
#' @param alpha a vector of size \code{2}.
#' @param x a matrix of size \code{c(nx, 2)} with angular coordinates.
#' @param N number of discretization points for approximating curve lengths.
#' Defaults to \code{5e2}.
#' @param L number of discretization points for computing the arc-length
#' parametrization curve lengths. Defaults to \code{5e2}.
#' @param der use derivatives to approximate curve lengths? Defaults to
#' \code{TRUE}.
#' @param shortest return the shortest possible distance? Defaults to
#' \code{TRUE}.
#' @param ridge_curve_grid if provided, the \code{ridge_curve} evaluated at
#' a grid of size \code{N}. If not provided, it is computed internally.
#' Useful for saving computations.
#' @param arclength use the arc-length parametrization to compute the
#' projections? This yields a more uniform grid for searching the projections.
#' Defaults to \code{TRUE}.
#' @inheritParams ridge_fourier_fit
#' @return
#' \itemize{
#'   \item \code{ridge_curve}: a matrix of size \code{c(nth, 2)} with the ridge
#'   curve evaluated at \code{theta}.
#'   \item \code{der_ridge_curve}: a matrix of size \code{c(nth, 2)} with the
#'   derivatives of the ridge curve evaluated at \code{theta}.
#'   \item \code{dist_ridge_curve}: the distance between two points along
#'   the ridge curve, a non-negative scalar.
#'   \item \code{proj_ridge_curve}: a list with (1) the \code{theta}'s that
#'   give the points in the ridge curve that are the closest (in the flat torus
#'   distance) to \code{x} (a matrix of size \code{c(nx, 2)}); (2) the indexes
#'   of \code{ridge_curve_grid} in which those \code{theta}'s were obtained.
#'   \item \code{arclength_ridge_curve}: a vector of size \code{N} giving the
#'   \code{theta} angles that yield a uniform-length grid of the ridge curve.
#' }
#' @examples
#' mu <- c(-0.5, 1.65)
#' th <- seq(-pi, pi, l = 200)
#' K <- 5
#' coefs <- list(cos_a = 1 / (1:(K + 1))^3, sin_b = 1 / (1:K)^3)
#' rid1 <- ridge_curve(theta = th, mu = mu, coefs = coefs, ind_var = 1)
#' rid2 <- ridge_curve(theta = th, mu = mu, coefs = coefs, ind_var = 2)
#' plot(mu[1], mu[2], xlim = c(-pi, pi), ylim = c(-pi, pi), axes = FALSE,
#'      xlab = expression(theta[1]), ylab = expression(theta[2]),
#'      pch = "*", col = 5, cex = 3)
#' sdetorus::linesTorus(rid1[, 1], rid1[, 2], col = 1)
#' sdetorus::linesTorus(rid2[, 1], rid2[, 2], col = 2)
#' abline(v = mu[1], lty = 3, col = 5)
#' abline(h = mu[2], lty = 3, col = 5)
#' points(ridge_curve(theta = mu[1], mu = mu, coefs = coefs, ind_var = 1),
#'        col = 1)
#' points(ridge_curve(theta = mu[2], mu = mu, coefs = coefs, ind_var = 2),
#'        col = 2)
#' sdetorus::torusAxis()
#'
#' ## der_ridge_curve
#'
#' th <- seq(-pi, pi, l = 10)
#' mu <- c(0.5, 1.5)
#' K <- 5
#' coefs <- list(cos_a = 1 / (1:(K + 1))^3, sin_b = 1 / (1:K)^3)
#' rid1 <- ridge_curve(theta = th, mu = mu, coefs = coefs, ind_var = 1)
#' rid2 <- ridge_curve(theta = th, mu = mu, coefs = coefs, ind_var = 2)
#' v1 <- der_ridge_curve(theta = th, mu = mu, coefs = coefs, ind_var = 1,
#'                       norm = 0.5)
#' v2 <- der_ridge_curve(theta = th, mu = mu, coefs = coefs, ind_var = 2,
#'                       norm = 0.5)
#' points(rid1, pch = 16, col = 1)
#' points(rid2, pch = 16, col = 2)
#' arrows(x0 = rid1[, 1], y0 = rid1[, 2],
#'        x1 = (rid1 + v1)[, 1], y1 = (rid1 + v1)[, 2],
#'        col = 3, angle = 5, length = 0.1)
#' arrows(x0 = rid2[, 1], y0 = rid2[, 2],
#'        x1 = (rid2 + v2)[, 1], y1 = (rid2 + v2)[, 2],
#'        col = 4, angle = 5, length = 0.1)
#'
#' ## dist_ridge_curve
#'
#' # Distances accuracy
#' a <- c(-pi / 2, pi)
#' mu <- c(-pi / 2, pi / 2)
#' dist_ridge_curve(alpha = a, mu = mu, coefs = coefs, der = TRUE, N = 1e6)
#' dist_ridge_curve(alpha = a, mu = mu, coefs = coefs, der = FALSE, N = 1e6)
#' dist_ridge_curve(alpha = a, mu = mu, coefs = coefs, der = TRUE, N = 1e2)
#' dist_ridge_curve(alpha = a, mu = mu, coefs = coefs, der = FALSE, N = 1e2)
#'
#' ## arclength_ridge_curve
#'
#' mu <- c(-pi / 2, pi / 2)
#' alpha <- arclength_ridge_curve(mu = mu, coefs = coefs, ind_var = 1, N = 25)
#' alpha <- sdetorus::toPiInt(c(alpha, alpha[1]))
#' rid <- ridge_curve(theta = alpha, mu = mu, coefs = coefs, ind_var = 1)
#' plot(mu[1], mu[2], pch = "*", col = 5, cex = 3, xlim = c(-pi, pi),
#'      ylim = c(-pi, pi), axes = FALSE, xlab = expression(theta[1]),
#'      ylab = expression(theta[2]))
#' sdetorus::linesTorus(rid[, 1], rid[, 2], col = 1, pch = 16)
#' points(rid[, 1], rid[, 2], pch = 16, col = 1)
#' abline(v = mu[1], lty = 3, col = 5)
#' abline(h = mu[2], lty = 3, col = 5)
#' sdetorus::torusAxis()
#'
#' ## proj_ridge_curve
#'
#' mu <- c(0, 0)
#' n <- 25
#' x <- matrix(runif(2 * n, -pi, pi), nrow = n, ncol = 2)
#' col <- rainbow(n)
#' th <- seq(-pi, pi, l = 100)
#' old_par <- par(no.readonly = TRUE)
#' par(mfrow = c(1, 2))
#' for (j in 1:2) {
#'
#'   plot(x, xlim = c(-pi, pi), ylim = c(-pi, pi), axes = FALSE,
#'        xlab = expression(theta[1]), ylab = expression(theta[2]), col = col)
#'   rid <- ridge_curve(theta = th, mu = mu, coefs = coefs, ind_var = j)
#'   sdetorus::linesTorus(x = rid[, 1], y = rid[, 2], lwd = 2)
#'   abline(v = mu[1], lty = 3)
#'   abline(h = mu[2], lty = 3)
#'   points(mu[1], mu[2], pch = "*", cex = 3)
#'   sdetorus::torusAxis()
#'   theta_projs <- proj_ridge_curve(x = x, mu = mu, coefs = coefs, ind_var = j,
#'                                   ridge_curve_grid = rid)$theta_proj
#'   projs <- ridge_curve(theta = theta_projs, mu = mu, coefs = coefs,
#'                        ind_var = j)
#'   points(projs, col = col, pch = 3)
#'   for (i in 1:n) {
#'
#'     sdetorus::linesTorus(x = c(x[i, 1], projs[i, 1]),
#'                          y = c(x[i, 2], projs[i, 2]), col = col[i], lty = 3)
#'
#'   }
#'
#' }
#' par(old_par)
#' @export
ridge_curve <- function(theta, mu = c(0, 0), coefs =
                          list(cos_a = c(0, 0), sin_b = 0),
                        ind_var = 1, at2 = TRUE) {

  # Independent variable can only be 1 or 2
  stopifnot(ind_var %in% 1:2)

  # Fourier atan2 or sine fit?
  if (at2) {

    # Check the length of coefficients
    stopifnot(is.list(coefs))
    stopifnot(length(coefs$cos_a) == (length(coefs$sin_b) + 1))

    # Fourier expansion for sine and cosine, using implicit column recycling
    k <- seq_along(coefs$sin_b)
    theta_cen <- c(0, theta - mu[ind_var])
    cos_k <- cos(k %o% theta_cen)
    sin_k <- sin(k %o% theta_cen)
    cos_fou <- coefs$cos_a[1] / 2 + colSums(cos_k * coefs$cos_a[-1])
    sin_fou <- colSums(sin_k * coefs$sin_b)
    r <- atan2(sin_fou, cos_fou)

    # Shift so r is centered at mu[ind_var]
    r <- sdetorus::toPiInt(mu[ifelse(ind_var == 1, 2, 1)] + r[-1] - r[1])

    # Return ridge curve depending on which variable represents theta
    return(unname(switch(ind_var, cbind(theta, r), cbind(r, theta))))

  } else {

    # Check coefficients
    stopifnot(is.numeric(coefs))

    # Fourier expansion, using implicit column recycling
    k <- seq_along(coefs)
    r <- colSums(sin(k %o% (theta - mu[ind_var])) * coefs)

    # Shift so r is centered at mu[ind_var]
    r <- sdetorus::toPiInt(mu[ifelse(ind_var == 1, 2, 1)] + r)

    # Return ridge curve depending on which variable represents theta
    return(unname(switch(ind_var, cbind(theta, r), cbind(r, theta))))

  }

}


#' @rdname ridge_curve
#' @export
der_ridge_curve <- function(theta, mu = c(0, 0), coefs =
                              list(cos_a = c(0, 0), sin_b = 0),
                            ind_var = 1, norm = NULL, at2 = TRUE) {

  # Independent variable can only be 1 or 2
  stopifnot(ind_var %in% 1:2)

  # Fourier atan2 or sine fit?
  if (at2) {

    # Check the length of coefficients
    stopifnot(is.list(coefs))
    stopifnot(length(coefs$cos_a) == (length(coefs$sin_b) + 1))

    # Fourier expansion for sine and cosine, using implicit column recycling
    k <- seq_along(coefs$sin_b)
    theta_cen <- theta - mu[ind_var]
    cos_k <- cos(k %o% theta_cen)
    sin_k <- sin(k %o% theta_cen)
    cos_fou <- coefs$cos_a[1] / 2 + colSums(cos_k * coefs$cos_a[-1])
    sin_fou <- colSums(sin_k * coefs$sin_b)

    # Derivatives
    d_cos_fou <- colSums(-sin_k * k * coefs$cos_a[-1])
    d_sin_fou <- colSums(cos_k * k * coefs$sin_b)
    d <- (-sin_fou * d_cos_fou + cos_fou * d_sin_fou) / (cos_fou^2 + sin_fou^2)

    # Derivative vectors depending on which variable represents theta
    v <- unname(switch(ind_var, cbind(1, d), cbind(d, 1)))

  } else {

    # Check coefficients
    stopifnot(is.numeric(coefs))

    # Derivative of the Fourier expansion, using implicit column recycling
    # Multiply by k due to chain rule. No shift required!
    k <- seq_along(coefs)
    d <- colSums(cos(k %o% (theta - mu[ind_var])) * k * coefs)

    # Derivative vectors depending on which variable represents theta
    v <- unname(switch(ind_var, cbind(1, d), cbind(d, 1)))

  }

  # Normalize vectors to a specific norm?
  if (!is.null(norm)) {

    v <- norm * v / sqrt(rowSums(v^2))

  }
  return(v)

}


#' @rdname ridge_curve
#' @export
dist_ridge_curve <- function(alpha, mu = c(0, 0), coefs =
                               list(cos_a = c(0, 0), sin_b = 0),
                             ind_var = 1, N = 5e2, der = TRUE,
                             shortest = TRUE, at2 = TRUE) {

  # Shortest sequence of alphas from alpha[1] to alpha[2]
  stopifnot(length(alpha) == 2)
  if (shortest) {

    alpha <- sdetorus::unwrapCircSeries(alpha)

  }
  alpha_seq <- seq(alpha[1], alpha[2], length.out = N)

  # Compute distance
  if (der) {

    # Approximate length by approximating the integral of the derivative norm
    der_curve <- der_ridge_curve(theta = alpha_seq, mu = mu, coefs = coefs,
                                 ind_var = ind_var, at2 = at2)
    d <- sdetorus::integrateSimp1D(fx = sqrt(rowSums(der_curve^2)),
                                   lengthInterval = abs(diff(alpha)))

  } else {

    # Approximate length by simply summing interpolating segments of the curve
    curve <- ridge_curve(theta = alpha_seq, mu = mu, coefs = coefs,
                         ind_var = ind_var, at2 = at2)
    d <- sum(sqrt(rowSums(sdetorus::diffCirc(curve)^2)))

  }
  return(d)

}


#' @rdname ridge_curve
#' @export
arclength_ridge_curve <- function(mu = c(0, 0), coefs =
                                    list(cos_a = c(0, 0), sin_b = 0),
                                  ind_var = 1, N = 5e2, L = 5e2, at2 = TRUE) {

  # Start measuring the distances from mu[ind_var]
  alpha <- mu[ind_var] +  seq(0, 2 * pi, l = L + 1)

  # Measure the lengths of the curve
  stopifnot(L >= N)
  alphas <- cbind(alpha[-(L + 1)], alpha[-1])
  alpha_distances <- apply(alphas, 1, function(a) {
     dist_ridge_curve(alpha = a, mu = mu, coefs = coefs, ind_var = ind_var,
                      N = 5e2, shortest = FALSE, der = TRUE, at2 = at2)
    })
  alpha_distances <- c(0, cumsum(alpha_distances))

  # Length of the curve
  max_length <- alpha_distances[L]

  # Uniform grid in the resulting lengths
  len_grid <- seq(0, max_length, l = N)

  # Interpolate to find the arcs that best give the lengths
  i <- findInterval(x = len_grid, vec = alpha_distances)
  w <- sdetorus::weightsLinearInterp1D(x = len_grid, g1 = alpha_distances[i],
                                       g2 = alpha_distances[i + 1])
  arcalpha <- rowSums(w * cbind(alpha[i], alpha[i + 1]))
  return(arcalpha)

}


#' @rdname ridge_curve
#' @export
proj_ridge_curve <- function(x, mu = c(0, 0), coefs =
                               list(cos_a = c(0, 0), sin_b = 0),
                             ind_var = 1, N = 5e2, ridge_curve_grid = NULL,
                             arclength = FALSE, at2 = TRUE) {

  # Compute ridge_curve_grid?
  if (is.null(ridge_curve_grid)) {

    if (arclength) {

      # Uniform grid in the resulting lengths
      theta <- arclength_ridge_curve(mu = mu, coefs = coefs, ind_var = ind_var,
                                     N = N, at2 = at2)

    } else {

      # Uniform grid in the parameters
      theta <- seq(-pi, pi, l = N + 1)[-(N + 1)]

    }

    # Create a grid
    ridge_curve_grid <- ridge_curve(theta = theta, mu = mu, coefs = coefs,
                                    ind_var = ind_var, at2 = at2)

  }

  # Closest points in ridge_curve_grid to x
  ind_grid <- apply(x, 1, function(th) {
    which.min(torus_dist(x = ridge_curve_grid, y = th, squared = TRUE))
  })
  theta_proj <- ridge_curve_grid[ind_grid, ind_var]
  return(list("theta_proj" = theta_proj, "ind_grid" = ind_grid))

}
