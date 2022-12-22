
#' @title Toroidal pairs plot
#'
#' @description Pairs plots for data on \eqn{[-\pi, \pi)^d}, \eqn{d\geq 2}.
#' The diagonal panels contain kernel density estimates tailored to
#' circular data.
#'
#' @inheritParams torus_dist
#' @param max_dim the maximum number of scores to produce the scores plot.
#' Defaults to \code{10}.
#' @param columns if specified, the variables to be plotted. If \code{NULL}
#' (the default), the first \code{max_dim} variables are plotted.
#' @param col_data color(s) for the data points. Defaults to \code{1}.
#' @param ylim_dens common \code{ylim} for the diagonal plots. Defaults to
#' \code{c(0, 1)}.
#' @param bwd type of bandwidth selector used in the kernel density plots.
#' Either \code{"ROT"}, \code{"EMI"}, \code{"AMI"}, \code{"LSCV"}, or
#' \code{"LCV"}. See \code{\link[DirStats:bw_dir_pi]{bw_dir_pi}} and
#' \code{\link[DirStats:bw_dir_pi]{bw_dir_cv}}. Defaults to \code{"ROT"}.
#' @param scales scales of the torus. Defaults to \code{rep(pi, ncol(x))}.
#' @details
#' The default bandwidth selector is the Rule-Of-Thumb (ROT) selector in
#' García-Portugués (2013). It is fast, yet it may oversmooth non-unimodal
#' densities. The EMI selector gives more flexible fits.
#' @return A \code{\link[ggplot2:ggplot2]{ggplot}}. The density plots show the
#' \link[=frechet_mean]{Fréchet means} (red bars) and the Fréchet standard
#' deviations (gray text).
#' @references
#' García-Portugués, E. (2013). Exact risk improvement of bandwidth selectors
#' for kernel density estimation with directional data. \emph{Electronic Journal
#' of Statistics}, 7:1655--1685. \doi{10.1214/13-ejs821}
#' @examples
#' # Generate data
#' n <- 50
#' set.seed(123456)
#' x <- sdetorus::toPiInt(rbind(
#'   mvtnorm::rmvnorm(n = n, mean = c(-pi, -pi) / 2,
#'                    sigma = diag(0.1, nrow = 2)),
#'   mvtnorm::rmvnorm(n = n, mean = c(-3 * pi / 2, 0) / 2,
#'           sigma = diag(0.1, nrow = 2)),
#'   mvtnorm::rmvnorm(n = n, mean = c(0, pi / 2),
#'                    sigma = diag(0.1, nrow = 2))
#' ))
#' col <- rainbow(3)[rep(1:3, each = n)]
#'
#' # Torus pairs
#' torus_pairs(x, col_data = col)
#' \donttest{
#' fit <- ridge_pca(x = x)
#' torus_pairs(fit$scores, col_data = col)}
#' @export
torus_pairs <- function(x, max_dim = 10, columns = NULL, col_data = 1,
                        ylim_dens = c(0, 1.5), bwd = "ROT",
                        scales = rep(pi, ncol(x))) {

  #Necessary variables
  d <- ncol(x)
  x <- as.data.frame(x)
  nam <- names(x)
  if (is.null(columns)) {

    columns <- 1:min(d, max_dim)

  }

  # Diagonal function
  kde_circle_pairs <- function(data, mapping, ...) {

    # Data index
    i <- match(x = as.character(mapping$x[[2]]), table = nam)

    # Half-period
    l <- scales[i]

    # Set labels for axes
    lims_l <- c(-l, l)
    brs_l <- seq(-l, l, length.out = 5)
    if (l == pi) {

      expr_l <- expression(-pi, -pi / 2, 0, pi / 2, pi)

    } else {

      expr_l <- expression(-L, -L / 2, 0, L / 2, L)

    }

    # Compute kde
    nth <- 500
    th <- seq(-pi, pi, l = nth)
    rad_data <- data[, i]
    cir_data <- DirStats::to_cir(rad_data / l * pi)
    h <- switch(bwd,
      "ROT" = DirStats::bw_dir_rot(data = cir_data),
      "AMI" = DirStats::bw_dir_ami(data = cir_data),
      "EMI" = DirStats::bw_dir_emi(data = cir_data)$h_opt,
      "LSCV" = DirStats::bw_dir_lscv(data = cir_data)$h_opt,
      "LCV" = DirStats::bw_dir_lcv(data = cir_data)$h_opt,
      stop("bwd must be \"ROT\", \"EMI\", \"AMI\", \"LSCV\", or \"LCV\".")
    )
    kde <- DirStats::kde_dir(x = DirStats::to_cir(th), data = cir_data, h = h)
    df <- data.frame("th" = th, "kde" = kde)

    # Plot
    p <- ggplot2::ggplot() +
      ggplot2::geom_line(mapping = ggplot2::aes(x = df$th * l / pi,
                                                y = df$kde * l / pi)) +
      ggplot2::scale_x_continuous(limits = lims_l, breaks = brs_l,
                                  labels = expr_l) +
      ggplot2::scale_y_continuous(limits = ylim_dens)

    # Variability
    frechet <- frechet_mean(x = rad_data, l = l)
    p <- p +
      ggplot2::geom_rug(mapping = ggplot2::aes(x = rad_data), alpha = 0.1,
                        size = 0.5) +
      ggplot2::geom_vline(xintercept = frechet$mu, colour = 2, alpha = 0.5) +
      ggplot2::annotate("text", x = 0, y = 0.6,
                        label = sprintf("%.2f", frechet$sd), size = 15,
                        color = gray(0.5, alpha = 0.5))
    return(p)

  }

  # Off-diagonal function
  points_torus_pairs <- function(data, mapping, ...) {

    # Data indexes
    ij <- match(x = as.character(c(mapping$x[[2]], mapping$y[[2]])),
                table = nam)

    # Half-periods
    li <- scales[ij[1]]
    lj <- scales[ij[2]]

    # Set labels for axes
    lims_i <- c(-li, li)
    lims_j <- c(-lj, lj)
    brs_i <- seq(-li, li, length.out = 5)
    brs_j <- seq(-lj, lj, length.out = 5)
    if (li == pi) {

      expr_i <- expression(-pi, -pi / 2, 0, pi / 2, pi)

    } else {

      expr_i <- expression(-L[i], -L[i] / 2, 0, L[i] / 2, L[i])

    }
    if (lj == pi) {

      expr_j <- expression(-pi, -pi / 2, 0, pi / 2, pi)

    } else {

      expr_j <- expression(-L[j], -L[j] / 2, 0, L[j] / 2, L[j])

    }

    # Regular plot
    p <- ggplot2::ggplot(data = data, mapping = mapping) +
      ggplot2::geom_point(..., colour = col_data, size = 0.5) +
      ggplot2::scale_x_continuous(limits = lims_i, breaks = brs_i,
                                  labels = expr_i) +
      ggplot2::scale_y_continuous(limits = lims_j, breaks = brs_j,
                                  labels = expr_j)
    return(p)

  }

  # Pairs plot
  p <- GGally::ggpairs(data = x[, columns],
                       lower = list(continuous = points_torus_pairs),
                       upper = list(continuous = points_torus_pairs),
                       diag = list(continuous = kde_circle_pairs),
                       columnLabels = nam)
  return(p)

}


#' @title Illustration of toroidal PCA via density ridges
#'
#' @description Shows the scores computation for PCA via density ridges on
#' \eqn{[-\pi, \pi)^2}.
#'
#' @param fit the output of \code{\link{ridge_pca}}.
#' @param n_max maximum number of data points to draw. These are sampled from
#' the data provided. Defaults to \code{500}.
#' @param projs draw projections? Defaults to \code{TRUE}.
#' @param projs_lines draw projection lines? Defaults to \code{TRUE}.
#' @param signs plot the original data points with \code{+} and \code{-} facets
#' depending on the signs of the second scores? Defaults to \code{TRUE}.
#' @inheritParams torus_pairs
#' @param col_projs a vector of size \code{2} giving the colors for the
#' curve-projected data and the ridge curve, respectively. Defaults to
#' \code{c(3, 4)}.
#' @param main caption of the plot. Empty by default.
#' @inheritParams ridge_curve
#' @inheritParams ridge_fourier_fit
#' @return Nothing, the functions are called to produce plots.
#' @examples
#' \donttest{
#' # Generate data
#' set.seed(987654321)
#' n <- 50
#' S1 <- rbind(c(1, -0.7), c(-0.7, 1))
#' S2 <- rbind(c(1, 0.5), c(0.5, 1))
#' x <- rbind(mvtnorm::rmvnorm(n, mean = c(0, pi / 2), sigma = S1),
#'            mvtnorm::rmvnorm(n, mean = c(pi, -pi / 2), sigma = S2))
#' x <- sdetorus::toPiInt(x)
#' col <- rainbow(2)[rep(1:2, each = n)]
#'
#' # ridge_pca and its visualization
#' fit <- ridge_pca(x = x, at2 = FALSE)
#' show_ridge_pca(fit = fit, col_data = col, at2 = FALSE)
#' fit2 <- ridge_pca(x = x, at2 = TRUE)
#' show_ridge_pca(fit = fit2, col_data = col, at2 = TRUE)}
#' @export
show_ridge_pca <- function(fit, n_max = 500, projs = TRUE, projs_lines = TRUE,
                           signs = TRUE, col_data = 1, col_projs = c(3, 4),
                           main = "", N = 5e2, at2 = TRUE) {

  # Extract and subsample data, if required
  x <- fit$data
  n <- nrow(x)
  ind <- sample(x = n, size = min(n_max, n))
  x <- x[ind, ]
  if (length(col_data) > 1) {

    col_data <- col_data[ind]

  }

  # Extract and scale the scores
  scores <- t(t(fit$scores[ind, ]) / fit$scales * pi)

  # Signs
  if (signs) {

    s <- sign(scores[, 2])
    pch <- c("-", "+")[(s + 1) / 2 + 1]

  } else {

    pch <- 16

  }

  # Plot data
  plot(x, xlim = c(-pi, pi), ylim = c(-pi, pi), axes = FALSE, col = col_data,
       xlab = expression(theta[1]), ylab = expression(theta[2]), pch = pch)
  title(main = main,
        sub = substitute("SS1 / SST = " * s,
                         list(s = sprintf("%.2f", fit$var_exp[1]))))
  sdetorus::torusAxis()
  abline(v = c(-pi, pi), col = "gray")
  abline(h = c(-pi, pi), col = "gray")

  # Plot ridge curve
  th_grid <- arclength_ridge_curve(mu = fit$mu_hat, coefs = fit$coefs_hat,
                                   ind_var = fit$ind_var, N = N, at2 = at2)
  th_grid <- c(th_grid, th_grid[1])
  y <- ridge_curve(theta = th_grid, mu = fit$mu_hat, coefs = fit$coefs_hat,
                   ind_var = fit$ind_var, at2 = at2)
  sdetorus::linesTorus(x = y[, 1], y = y[, 2], col = col_projs[2], lwd = 2)

  # Projections and lines
  if (projs || projs_lines) {

    ridge_projs <- ridge_curve(theta = scores[, 1] +
                                 fit$mu_hat[fit$ind_var],
                               mu = fit$mu_hat, coefs = fit$coefs_hat,
                               ind_var = fit$ind_var, at2 = at2)
    if (projs) {

      points(ridge_projs, cex = 0.5, pch = 16, col = col_projs[1])

    }
    if (projs_lines) {

      for (i in seq_along(ind)) {

        sdetorus::linesTorus(x = c(x[i, 1], ridge_projs[i, 1]),
                             y = c(x[i, 2], ridge_projs[i, 2]),
                             col = ggplot2::alpha(colour = col_projs[1],
                                                  alpha = 0.25))

      }

    }

  }

  # Backwards mean and antipodal backwards mean
  points(fit$mu_hat[1], fit$mu_hat[2], pch = "*", cex = 4, col = col_projs[2])

}
