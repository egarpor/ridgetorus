
#' @title Toroidal PCA via density ridges
#'
#' @description This function computes the whole process of toroidal PCA
#' via density ridges on a given sample: parameter estimation of the
#' underlying distribution, estimation of the connected component of the ridge,
#' and determination of its Fourier expansion from which to obtain the first
#' and second scores.
#'
#' @param x matrix of dimension  \code{c(n, 2)} containing the \code{n}
#' observations of the pair of angles.
#' @param type either \code{"bvm"} (bivariate sine von Mises), \code{"bwc"}
#' (bivariate wrapped Cauchy), or \code{"auto"} (default). \code{"auto"}
#' performs both fits and uses the one with lowest BIC.
#' @inheritParams ridge_curve
#' @inheritParams ridge_scores
#' @param lrts run \code{\link{biv_lrt}} to check the null hypothesis of
#' homogeneous concentration parameters using likelihood ratio tests? If
#' \code{TRUE} (default), enforces the special horizontal/vertical/diagonal
#' cases to become "sticky" fits.
#' @param alpha significance level for the homogeneity test.
#' @inheritParams biv_lrt
#' @inheritParams ridge_fourier_fit
#' @return A list with:
#' \item{mu_hat}{estimated circular means of the sample.}
#' \item{coefs_hat}{estimated Fourier coefficients.}
#' \item{ind_var}{indexing variable.}
#' \item{scores}{scores for each of the sample points.}
#' \item{var_exp}{percentage of explained variance.}
#' \item{fit_mle}{maximum likelihood fit.}
#' \item{bic_fit}{BIC of the fit.}
#' \item{data}{original sample.}
#' \item{scales}{vector of length 2 with the scale limits for the axes.}
#' \item{type}{type of fit performed.}
#' \item{p_hom}{\eqn{p}-value of the homogeneity test.}
#' \item{p_indep}{\eqn{p}-value of the independence test.}
#' @examples
#' \donttest{
#' ## Bivariate von Mises
#'
#' n <- 100
#' x <- r_bvm(n = n, mu = c(1, 2), kappa = c(0.4, 0.4, 0.5))
#' fit <- ridge_pca(x = x, type = "bvm")
#' show_ridge_pca(fit = fit, col_data = "red")
#'
#' x <- r_bvm(n = n, mu = c(2, 1), kappa = c(1, 2, 0))
#' fit <- ridge_pca(x = x, type = "bvm")
#' show_ridge_pca(fit = fit, col_data = "red")
#'
#' x <- r_bvm(n = n, mu = c(2, 1), kappa = c(3, 2, 0))
#' fit <- ridge_pca(x = x, type = "bvm")
#' show_ridge_pca(fit = fit, col_data = "red")
#'
#' ## Bivariate wrapped Cauchy
#'
#' x <- r_bwc(n = n, mu = c(1, 2), xi = c(0.2, 0.2, 0.5))
#' fit <- ridge_pca(x = x, type = "bwc")
#' show_ridge_pca(fit = fit, col_data = "red")
#'
#' x <- r_bwc(n = n, mu = c(1, 2), xi = c(0.2, 0.8, 0))
#' fit <- ridge_pca(x = x, type = "bwc")
#' show_ridge_pca(fit = fit, col_data = "red")
#'
#' x <- r_bwc(n = n, mu = c(1, 2), xi = c(0.5, 0.2, 0))
#' fit <- ridge_pca(x = x, type = "bwc")
#' show_ridge_pca(fit = fit, col_data = "red")}
#' @export
ridge_pca <- function(x, type = c("auto", "bvm", "bwc")[1], N = 5e2, K = 15,
                      scale = TRUE, lrts = TRUE, alpha = 0.05, at2 = TRUE,
                      ...) {

  # Initial checks
  if (!(is.matrix(x) || is.data.frame(x))) {

    stop("x must be a two column matrix or dataframe.")

  }
  if (ncol(x) != 2) {

    stop("x must be a two column matrix or dataframe.")

  }
  if (!(type %in% c("auto", "bvm", "bwc"))) {

    stop("Type is one of: auto, bvm or bwc.")

  }
  if (N <= 0) {

    stop(paste("Introduce a positive number N of discretization points",
               "for approximating curve lengths."))

  }
  if (K <= 0) {

    stop("Introduce a positive number K of Fourier expansion terms.")

  }
  if (alpha <= 0 || alpha > 1) {

    stop("The significance level alpha must be between 0 and 1.")

  }
  if (!is.logical(scale)) {

    stop("scale must be a logical value: TRUE or FALSE.")

  }
  if (!is.logical(lrts)) {

    stop("hom must be a logical value: TRUE or FALSE.")

  }
  if (!is.logical(at2)) {

    stop("at2 must be a logical value: TRUE or FALSE.")

  }

  # Remove missing data
  if (anyNA(x)) {

    message("Removing missing data")
    x <- na.omit(x)

  }

  # Sample size
  n <- nrow(x)

  # Estimation of underlying distribution
  if (type == "auto") {

    fit_bvm <- fit_bvm_mle(x = x, ...)
    fit_bwc <- fit_bwc_mle(x = x, ...)
    bvm_bwc <- (fit_bwc$opt$value < fit_bvm$opt$value) + 1
    fit_mle <- switch(bvm_bwc, fit_bvm, fit_bwc)
    type <- switch(bvm_bwc, "bvm", "bwc")

  } else {

    fit_mle <- switch(type,
                      "bvm" = fit_bvm_mle(x = x, ...),
                      "bwc" = fit_bwc_mle(x = x, ...),
                      stop("type must be \"auto\", \"bvm\" or \"bwc\"."))

  }

  # Likelihood ratio tests to make the special horizontal/vertical/diagonal
  # ridges to be "sticky" fits
  pval_indep <- pval_hom <- NA
  fit_mle0 <- fit_mle
  if (lrts) {

    # Test independence
    lrt_indep <- biv_lrt(x = x, fit_mle = fit_mle0, type = type, indep = TRUE,
                         hom = FALSE, ...)
    pval_indep <- lrt_indep$p.value
    indep <- pval_indep >= alpha
    if (indep) {

      fit_mle0 <- lrt_indep$fit_null

    }

    # Test homogeneity
    lrt_hom <- biv_lrt(x = x, fit_mle = fit_mle0, type = type, indep = indep,
                       hom = TRUE, ...)
    pval_hom <- lrt_hom$p.value
    hom <- pval_hom >= alpha
    if (hom) {

      fit_mle0 <- lrt_hom$fit_null

    }
    df <- lrt_hom$df

    # Undefined ridge?
    if (hom && indep) {

      warning(paste("The ridge is not well-defined: homogeneity and",
                    "independence is not rejected, so the fit is unstable.",
                    "Skipping horizontal/vertical/diagonal ridge."))
      fit_mle0 <- fit_mle
      df <- 0

    }

  }

  # Common values
  bic_fit <- ifelse(lrts, 5 - df, 5) * log(n) + 2 * fit_mle0$opt$value
  mu_hat <- c(fit_mle0$mu1, fit_mle0$mu2)
  conc_hat <- unname(unlist(fit_mle0[3:5]))

  # Connected component of the ridge
  if (type == "bvm") {

    ridge <- ridge_bvm(mu = mu_hat, kappa = conc_hat,
                       subint_1 = 5e2, subint_2 = 5e2)

  } else if (type == "bwc") {

    ridge <- ridge_bwc(mu = mu_hat, xi = conc_hat,
                       subint_1 = 5e2, subint_2 = 5e2)

  }

  # Fourier fit indexed by lowest kappa or xi variable
  ind_var <- order(conc_hat[1:2])[1]
  if (ind_var == 2) {

    ridge <- sdetorus::toPiInt(t(t(ridge[, 2:1]) - mu_hat[2:1]))

  } else {

    ridge <- sdetorus::toPiInt(t(t(ridge) - mu_hat))

  }

  # Estimate coefficients
  coefs_hat <- ridge_fourier_fit(curve = ridge, K = K, at2 = at2)

  # Scores computation
  scores <- ridge_scores(x = x, mu = mu_hat, coefs = coefs_hat,
                         ind_var = ind_var, N = N, scale = scale, at2 = at2)

  # Proportion of explained variances
  var_exp <- frechet_ss(x = scores$scores, l = scores$scales)$var_exp

  # Ridge fit
  return(list("mu_hat" = mu_hat, "coefs_hat" = coefs_hat, "ind_var" = ind_var,
              "scores" = scores$scores, "var_exp" = var_exp,
              "fit_mle" = fit_mle0, "bic_fit" = bic_fit, "data" = x,
              "scales" = scores$scales, "type" = type, "p_hom" = pval_hom,
              "p_indep" = pval_indep))

}
