
#' @title Scores and scales for Fourier-fitted ridge curves
#'
#' @description Computation of PCA scores for \link[=ridge_curve]{
#' Fourier-fitted ridge curves}. The scores are defined as follows:
#' \itemize{
#'   \item First scores: signed distances along the ridge curve of the data
#'   projections to \eqn{\mu}.
#'   \item Second scores: signed toroidal distances from the data points to
#'   their ridge projections.
#' }
#' The scores can be scaled to \eqn{(-\pi, \pi)} or remain as
#' \eqn{(l / 2, m_2)}, where \eqn{l} is the length of the curve and \eqn{m_2}
#' is the maximal absolute second score.
#'
#' @inheritParams ridge_curve
#' @param scale scale the resulting scores to \eqn{[-\pi, \pi)^2}? Defaults
#' to \code{TRUE}.
#' @inheritParams ridge_fourier_fit
#' @param L grid along he variable \code{ind_var} used for searching the
#' maximum allowed second score. Defaults to \code{25}.
#' @param f factor for shrinking the grid on the variable that is different to
#' \code{ind_var}. Defaults to \code{2}.
#' @details
#' The mean \eqn{\mu} corresponds to the first score being null.
#' @return \code{ridge_scores} returns a list with:
#' \item{scores}{a matrix of size \code{c(nx, 2)} with the ridge scores.}
#' \item{scales}{a vector of length 2 with the scale limits for the axes.}
#' \code{max_score_2} computes the maximum allowed second score to rescale if
#' \code{scale = TRUE}.
#' @examples
#' mu <- c(-0.5, 1.65)
#' th <- seq(-pi, pi, l = 200)
#' K <- 5
#' coefs <- list(cos_a = 1 / (1:(K + 1))^3, sin_b = 1 / (1:K)^3)
#' n <- 10
#' col <- rainbow(n)
#'
#' set.seed(13213)
#' old_par <- par(no.readonly = TRUE)
#' par(mfrow = c(2, 2))
#' for (j in 1:2) {
#'
#'   # Simulate synthetic data close to the ridge curve
#'   rid <- ridge_curve(theta = th, mu = mu, coefs = coefs, ind_var = j)
#'   ind <- sort(sample(length(th), size = n))
#'   eps <- 0.25 * matrix(runif(2 * n, -1, 1), nrow = n, ncol = 2)
#'   x <- sdetorus::toPiInt(rid[ind, ] + eps)
#'
#'   # Plot ridge and synthetic data, with signs from the second scores
#'   s <- ridge_scores(x, mu = mu, coefs = coefs, ind_var = j)$scores
#'   plot(x, xlim = c(-pi, pi), ylim = c(-pi, pi), axes = FALSE,
#'        xlab = expression(theta[1]), ylab = expression(theta[2]), col = col,
#'        pch = ifelse(sign(s[, 2]) == 1, "+", "-"), cex = 1.25)
#'   sdetorus::linesTorus(rid[, 1], rid[, 2], lwd = 2)
#'   abline(v = mu[1], lty = 3)
#'   abline(h = mu[2], lty = 3)
#'   points(mu[1], mu[2], pch = "*", cex = 3)
#'   sdetorus::torusAxis()
#'
#'   # Projections
#'   theta_projs <- proj_ridge_curve(x = x, mu = mu, coefs = coefs,
#'                                   ind_var = j, ridge_curve_grid = rid,
#'                                   )$theta_proj
#'   projs <- ridge_curve(theta = theta_projs, mu = mu, coefs = coefs,
#'                        ind_var = j)
#'   for (i in 1:n) {
#'
#'     sdetorus::linesTorus(x = c(x[i, 1], projs[i, 1]),
#'                          y = c(x[i, 2], projs[i, 2]),
#'                          col = col[i], lty = 3)
#'
#'   }
#'
#'   # Scores plot
#'   plot(s, xlim = c(-pi, pi), ylim = c(-pi, pi), axes = FALSE,
#'        xlab = "Score 1", ylab = "Score 2", col = col,
#'        pch = ifelse(sign(s[, 2]) == 1, "+", "-"))
#'   sdetorus::torusAxis()
#'   abline(v = 0, lty = 3)
#'   abline(h = 0, lty = 3)
#'   points(0, 0, pch = "*", cex = 3)
#'
#' }
#' par(old_par)
#' @export
ridge_scores <- function(x, mu = c(0, 0), coefs =
                           list(cos_a = c(0, 0), sin_b = 0),
                         ind_var = 1, N = 5e2, scale = TRUE, at2 = TRUE) {

  ## First scores: distances, along the ridge curve, of the data projections
  ## to mu

  # Construct arc-length parametrization and ridge_curve_grid
  alpha <- arclength_ridge_curve(mu = mu, coefs = coefs, ind_var = ind_var,
                                 L = 10 * N, N = N, at2 = at2)
  ridge_curve_grid <- ridge_curve(theta = alpha, mu = mu, coefs = coefs,
                                  ind_var = ind_var, at2 = at2)

  # Projections, including mu for later recentering (so score_1 = 0 if and
  # only if x = mu)
  projs_curve <- proj_ridge_curve(x = rbind(mu, x), mu = mu, coefs = coefs,
                                  ind_var = ind_var, N = N, ridge_curve_grid =
                                    ridge_curve_grid, at2 = at2)
  scores_1 <- alpha[projs_curve$ind_grid]
  scores_1_mu <- scores_1[1]
  scores_1 <- scores_1[-1]

  ## Second scores: *signed* distances to ridge curve

  # Unsigned distances
  y <- ridge_curve(theta = scores_1, mu = mu, coefs = coefs, ind_var = ind_var,
                   at2 = at2)
  scores_2 <- torus_dist(x = x, y = y)

  # Assign signs using the sign of the angle between the tangent and normal
  # vectors at the projections
  norms <- sdetorus::toPiInt(y - x)
  tangents <- der_ridge_curve(theta = scores_1, mu = mu, coefs = coefs,
                              ind_var = ind_var, at2 = at2)
  signs_scores_2 <- sign(sdetorus::toPiInt(
    atan2(tangents[, 2], tangents[, 1]) - atan2(norms[, 2], norms[, 1])
  ))
  scores_2 <- signs_scores_2 * scores_2

  # Centering of scores_1
  scores_1 <- sdetorus::toPiInt(scores_1 - scores_1_mu)

  ## Scale the scores?

  max_sc2 <- max_score_2(mu = mu, coefs = coefs, ind_var = ind_var,
                         L = 25, f = 2, at2 = at2)
  if (scale) {

    # Scale second scores to lie in (-pi, pi) (analogously as it is done with
    # the scores_1 via the arc-length parametrization of the ridge curve)
    scores_2 <- scores_2 * (pi / max_sc2)
    scales <- c(pi, pi)

  } else {

    # Unscale first scores
    length_ridge <- dist_ridge_curve(alpha = c(0, 2 * pi), mu = mu,
                                     coefs = coefs, ind_var = ind_var,
                                     N = 1e4, shortest = FALSE, at2 = at2)
    scores_1 <- scores_1 * length_ridge / (2 * pi)
    scales <- c(length_ridge / 2, max_sc2)

  }

  # Return scores and scales
  return(list("scores" = cbind(scores_1, scores_2), "scales" = scales))

}


#' @rdname ridge_scores
#' @export
max_score_2 <- function(mu = c(0, 0), coefs =
                          list(cos_a = c(0, 0), sin_b = 0),
                        ind_var = 1, L = 25, f = 2, at2 = TRUE) {

  # Full ridge grid
  th_1 <- seq(-pi, pi, l = f * L + 1)[-(f * L + 1)]
  ridge_curve_grid <- ridge_curve(theta = th_1, mu = mu, coefs = coefs,
                                  ind_var = ind_var, at2 = at2)

  # Grid targeting the regions giving the maximal separation from ridge curve
  th_2 <- sdetorus::toPiInt(pi + seq(-pi / f, pi / f, l = L))
  sc_2_grid <- ridge_curve(theta = th_1, mu = mu, coefs = coefs,
                           ind_var = ind_var, at2 = at2)
  not_ind_var <- ifelse(ind_var == 1, 2, 1)
  sc_2_grid <- c(sdetorus::toPiInt(outer(sc_2_grid[, not_ind_var], th_2, "+")))
  sc_2_grid <- switch(ind_var, cbind(th_1, sc_2_grid), cbind(sc_2_grid, th_1))

  # Projections
  proj_scores_2_grid <- proj_ridge_curve(x = sc_2_grid,
                                         ridge_curve_grid = ridge_curve_grid,
                                         ind_var = ind_var)
  proj_scores_2_grid <- ridge_curve_grid[proj_scores_2_grid$ind_grid, ]

  # Maximum distance
  max_score_2 <- max(torus_dist(x = sc_2_grid, y = proj_scores_2_grid,
                                squared = FALSE))
  return(max_score_2)

}
