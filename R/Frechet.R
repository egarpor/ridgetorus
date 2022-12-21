
#' @title Fréchet statistics on the torus
#'
#' @description Computes the Fréchet mean, variance, and standard deviation of
#' a sample on the \eqn{d}-torus \eqn{[-l, l)^d}, \eqn{d\geq 1}, with
#' \eqn{-l \equiv l} identified.
#'
#' @param x sample of angles on \eqn{[-l, l)}, a vector or a matrix.
#' @param l half-period of the circular data. Can be a vector of length
#' \code{ncol(x)} if \code{x} is a matrix. Defaults to \code{pi}.
#' @param N size of the grid in \eqn{[-l, l)} for the exhaustive search of
#' the mean. Defaults to \code{5e2}.
#' @param draw_plot draw a diagnostic plot showing the Fréchet loss function?
#' Defaults to \code{FALSE}.
#' @return
#' \itemize{
#'   \item \code{frechet_mean}: a list with the marginal Fréchet means
#'    (\code{mu}), variances (\code{var}), and standard deviations (\code{sd}).
#'   \item \code{frechet_ss}: a list with the Fréchet variances (\code{var})
#'   and the cumulative proportion of total variance explained (\code{var_exp}).
#' }
#' @examples
#' ## Circular data
#'
#' # Sample from a wrapped normal
#' x <- sdetorus::toPiInt(rnorm(n = 5e2, mean = 2, sd = 1))
#' frechet_mean(x = x)
#'
#' # Sample from a bimodal distribution
#' x <- sdetorus::toPiInt(rnorm(n = 5e2, mean = c(1, -2), sd = c(0.5, 0.75)))
#' frechet_mean(x = x)
#'
#' # Periodic data in [-2, 2)
#' x <- sdetorus::toInt(rnorm(n = 5e2, mean = c(-2, 1), sd = 2:1),
#'                      a = -2, b = 2)
#' frechet_mean(x = x, l = 2)
#'
#' ## Toroidal data
#'
#' # Sample from a multivariate wrapped normal
#' n <- 50
#' S <- rbind(c(2.5, -0.2, 0.5),
#'            c(-0.2, 1.5, -0.5),
#'            c(0.5, -0.5, 0.75))
#' x <- sdetorus::toPiInt(mvtnorm::rmvnorm(n, mean = c(0, 1.5, -2), sigma = S))
#' (f <- frechet_mean(x = x))
#'
#' # Total Fréchet variance is sum of marginal variances
#' sum(torus_dist(x, y = f$mu, squared = TRUE)) / n
#' sum(f$var)
#'
#' # Cumulative proportion of variances
#' frechet_ss(x)
#' @name frechet
#' @rdname frechet
#' @export
frechet_mean <- function(x, l = pi, N = 5e2, draw_plot = FALSE) {

  # As matrix if x is a vector (circular data)
  if (is.vector(x)) {

    x <- matrix(x, ncol = 1)

  }

  # Create n, d, and l
  n <- nrow(x)
  d <- ncol(x)
  if (length(l) == 1) {

    l <- rep(l, d)

  }

  # Circular Fréchet mean, variance, and standard deviation
  frechet_cir <- function(j) {

    # Exhaustive search in a grid of size N
    th_grid <- seq(-l[j], l[j], length.out = N + 1)[- (N + 1)]
    ss <- sapply(th_grid, function(th) {
      sum(sdetorus::toInt(th - x[, j], a = -l[j], b = l[j])^2)
      }) / n

    # Mean and variance
    ind <- which.min(ss)
    mu <- th_grid[ind]
    var <- ss[ind]

    # Draw diagnostic plots?
    if (draw_plot) {

      # Minimum
      plot(th_grid, ss, type = "l", axes = FALSE, main = paste("Variable", j),
           xlab = substitute(x[jj], list(jj = j)),
           ylab = "Fr\u00e9chet variance")
      abline(v = mu, col = 2, lwd = 2)
      if (l[j] == pi) {

        sdetorus::torusAxis(1)

      } else {

        axis(1)
        box()

      }
      axis(2)

    }

    # Return Fréchet mean, variance, and standard deviation
    return(c(mu, var, sqrt(var)))

  }

  # Fréchet analysis
  frechet_torus <- sapply(1:d, frechet_cir)

  # Return Fréchet mean, variance, and standard deviation
  return(list("mu" = frechet_torus[1, ], "var" = frechet_torus[2, ],
              "sd" = frechet_torus[3, ]))

}


#' @rdname frechet
#' @export
frechet_ss <- function(x, l = pi, N = 5e2, draw_plot = FALSE) {

  vars <- frechet_mean(x = x, l = l, N = N, draw_plot = draw_plot)$var
  return(list("var" = vars, "var_exp" = cumsum(vars) / sum(vars)))

}


#' @title Toroidal distances
#'
#' @description Computation of distances on \eqn{[-\pi, \pi)^d}, \eqn{d\geq 1},
#' between two sets of observations.
#'
#' @param x a matrix of size \code{c(nx, d)} with angles on \eqn{[-\pi, \pi)}.
#' @param y either a matrix with the same size as \code{x} or a vector of
#' size \code{nx}.
#' @param squared return the squared distance? Defaults to \code{FALSE}.
#' @details
#' The maximal distance on \eqn{[-\pi, \pi)^d} is \eqn{\sqrt{d}\pi}.
#' @return A vector of size \code{nx} with the distances between the
#' observations of \code{x} and \code{y}.
#' @examples
#' # Illustration of torus distances
#' n <- 10
#' x <- matrix(runif(2 * n, -pi, pi), nrow = n, ncol = 2)
#' y <- c(pi / 2, pi / 3)
#' col <- rainbow(n)
#' plot(x, xlim = c(-pi, pi), ylim = c(-pi, pi), axes = FALSE, col = col,
#'      xlab = expression(theta[1]), ylab = expression(theta[2]), pch = 16)
#' sdetorus::torusAxis()
#' points(y[1], y[2], col = 1, pch = 17)
#' for (i in 1:n) {
#'
#'   sdetorus::linesTorus(x = c(x[i, 1], y[1]),
#'                        y = c(x[i, 2], y[2]), lty = 2, col = col[i])
#'
#' }
#' text(x = x, labels = sprintf("%.2f", torus_dist(x, y)), col = col, pos = 1)
#' @export
torus_dist <- function(x, y, squared = FALSE) {

  stopifnot(is.data.frame(x) || is.matrix(x))
  stopifnot(ncol(x) == length(y) || identical(dim(x), dim(y)))
  d2 <- switch(is.vector(y) + 1,
               rowSums(sdetorus::toPiInt(x - y)^2),
               colSums(sdetorus::toPiInt(t(x) - y)^2))
  return(switch(squared + 1, sqrt(d2), d2))

}
