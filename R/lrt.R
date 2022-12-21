
#' @title Tests of homogeneity and independence in bivariate sine von
#' Mises and wrapped Cauchy distributions
#'
#' @description Performs the following likelihood ratio tests for the
#' concentrations in bivariate sine von Mises and wrapped Cauchy distributions:
#' (1) \emph{homogeneity}: \eqn{H_0:\kappa_1=\kappa_2} vs.
#' \eqn{H_1:\kappa_1\neq\kappa_2}, and \eqn{H_0:\xi_1=\xi_2} vs.
#' \eqn{H_1:\xi_1\neq\xi_2}, respectively;
#' (2) \emph{independence}: \eqn{H_0:\lambda=0} vs.
#' \eqn{H_1:\lambda\neq0}, and \eqn{H_0:\rho=0} vs. \eqn{H_1:\rho\neq0}.
#' The tests (1) and (2) can be performed simultaneously.
#'
#' @inheritParams bwc
#' @param hom test the homogeneity hypothesis? Defaults to \code{FALSE}.
#' @param indep test the independence hypothesis? Defaults to \code{FALSE}.
#' @param fit_mle output of \code{\link{fit_bvm_mle}} or
#' \code{\link{fit_bwc_mle}} with \code{hom = FALSE}. Computed internally if
#' not provided.
#' @param type either \code{"bvm"} (bivariate sine von Mises) or \code{"bwc"}
#' (bivariate wrapped Cauchy).
#' @inheritParams ridge_pca
#' @param ... optional parameters passed to \code{\link{fit_bvm_mle}} and
#' \code{\link{fit_bwc_mle}}, such as \code{start}, \code{lower}, or
#' \code{upper}.
#' @return A list with class \code{htest}:
#' \item{statistic}{the value of the likelihood ratio test statistic.}
#' \item{p.value}{the \eqn{p}-value of the test (computed using the asymptotic
#' distribution).}
#' \item{alternative}{a character string describing the alternative
#' hypothesis.}
#' \item{method}{description of the type of test performed.}
#' \item{df}{degrees of freedom.}
#' \item{data.name}{a character string giving the name of \code{theta}.}
#' \item{fit_mle}{maximum likelihood fit.}
#' \item{fit_null}{maximum likelihood fit under the null hypothesis.}
#' @references
#' Kato, S. and Pewsey, A. (2015). A Möbius transformation-induced distribution
#' on the torus. \emph{Biometrika}, 102(2):359--370. \doi{10.1093/biomet/asv003}
#'
#' Singh, H., Hnizdo, V., and Demchuk, E. (2002). Probabilistic model for two
#' dependent circular variables. \emph{Biometrika}, 89(3):719--723.
#' \doi{10.1093/biomet/89.3.719}
#' @examples
#' ## Bivariate sine von Mises
#'
#' # Homogeneity
#' n <- 200
#' mu <- c(0, 0)
#' kappa_0 <- c(1, 1, 0.5)
#' kappa_1 <- c(0.7, 0.1, 0.25)
#' samp_0 <- r_bvm(n = n, mu = mu, kappa = kappa_0)
#' samp_1 <- r_bvm(n = n, mu = mu, kappa = kappa_1)
#' biv_lrt(x = samp_0, hom = TRUE, type = "bvm")
#' biv_lrt(x = samp_1, hom = TRUE, type = "bvm")
#'
#' # Independence
#' kappa_0 <- c(0, 1, 0)
#' kappa_1 <- c(1, 0, 1)
#' samp_0 <- r_bvm(n = n, mu = mu, kappa = kappa_0)
#' samp_1 <- r_bvm(n = n, mu = mu, kappa = kappa_1)
#' biv_lrt(x = samp_0, indep = TRUE, type = "bvm")
#' biv_lrt(x = samp_1, indep = TRUE, type = "bvm")
#'
#' # Independence and homogeneity
#' kappa_0 <- c(3, 3, 0)
#' kappa_1 <- c(3, 1, 0)
#' samp_0 <- r_bvm(n = n, mu = mu, kappa = kappa_0)
#' samp_1 <- r_bvm(n = n, mu = mu, kappa = kappa_1)
#' biv_lrt(x = samp_0, indep = TRUE, hom = TRUE, type = "bvm")
#' biv_lrt(x = samp_1, indep = TRUE, hom = TRUE, type = "bvm")
#'
#' ## Bivariate wrapped Cauchy
#'
#' # Homogeneity
#' xi_0 <- c(0.5, 0.5, 0.25)
#' xi_1 <- c(0.7, 0.1, 0.5)
#' samp_0 <- r_bwc(n = n, mu = mu, xi = xi_0)
#' samp_1 <- r_bwc(n = n, mu = mu, xi = xi_1)
#' biv_lrt(x = samp_0, hom = TRUE, type = "bwc")
#' biv_lrt(x = samp_1, hom = TRUE, type = "bwc")
#'
#' # Independence
#' xi_0 <- c(0.1, 0.5, 0)
#' xi_1 <- c(0.3, 0.5, 0.2)
#' samp_0 <- r_bwc(n = n, mu = mu, xi = xi_0)
#' samp_1 <- r_bwc(n = n, mu = mu, xi = xi_1)
#' biv_lrt(x = samp_0, indep = TRUE, type = "bwc")
#' biv_lrt(x = samp_1, indep = TRUE, type = "bwc")
#'
#' # Independence and homogeneity
#' xi_0 <- c(0.2, 0.2, 0)
#' xi_1 <- c(0.1, 0.2, 0.1)
#' samp_0 <- r_bwc(n = n, mu = mu, xi = xi_0)
#' samp_1 <- r_bwc(n = n, mu = mu, xi = xi_1)
#' biv_lrt(x = samp_0, indep = TRUE, hom = TRUE, type = "bwc")
#' biv_lrt(x = samp_1, indep = TRUE, hom = TRUE, type = "bwc")
#' @export
biv_lrt <- function(x, hom = FALSE, indep = FALSE, fit_mle = NULL, type, ...) {

  # Ensure a proper lrt is conducted
  if (!hom && !indep) {

    stop("At least one of hom or indep must be TRUE.")

  }

  # Fit
  fit_type_mle <- switch(type, "bvm" = fit_bvm_mle, "bwc" = fit_bwc_mle,
                         stop("type must be \"bvm\" or \"bwc\"."))

  # Fit unconstrained and constrained models
  if (is.null(fit_mle)) {

    fit_mle <- fit_type_mle(x = x, hom = FALSE, indep = FALSE, ...)

  }
  op <- list(...)
  op$start <- unlist(fit_mle[1:5])
  fit_null <- do.call(what = "fit_type_mle",
                      args = c(list(x = x, hom = hom, indep = indep), op))
  minus_loglik <- fit_mle$opt$value
  minus_loglik_null <- fit_null$opt$value

  # Is the constrained model better than the unconstrained? If so, improve
  # the unconstrained model
  if (minus_loglik_null < minus_loglik) {

    message(paste0(type, "_loglik_null(hom = ", hom, ", indep = ", indep,
                   ") > ", type, "_loglik, refitting unconstrained model"))
    start <- unname(unlist(fit_null[1:5]))
    fit_mle <- fit_type_mle(x = x, hom = FALSE, indep = FALSE,
                            start = start, ...)
    minus_loglik <- fit_mle$opt$value

  }

  # Likelihood ratio test
  lr_stat <- 2 * (minus_loglik_null - minus_loglik)
  df <- hom + indep
  pvalue <- pchisq(q = lr_stat, df = df, lower.tail = FALSE)

  # Test information
  method <- "Likelihood ratio test for H0:"
  if (hom) {

    if (indep) {

      dep <- ifelse(type == "bvm", "lambda", "rho")
      method <-
        paste(method, ifelse(type == "bvm", "kappa1 = kappa2 and",
                             "xi1 = xi2 and"), dep, "= 0 in a bivariate",
              ifelse(type == "bvm", "sine von Mises", "wrapped Cauchy"))
      alt <- paste(ifelse(type == "bvm", "H1: kappa1 != kappa2 or",
                          "H1: xi1 != xi2 or"), dep, "!= 0")

    } else {

      method <- paste(method, ifelse(type == "bvm", "kappa1 = kappa2",
                                     "xi1 = xi2"), "in a bivariate",
                      ifelse(type == "bvm", "sine von Mises", "wrapped Cauchy"))
      alt <- ifelse(type == "bvm", "H1: kappa1 != kappa2", "H1: xi1 != xi2")

    }

  } else {

    dep <- ifelse(type == "bvm", "lambda", "rho")
    method <- paste(method, dep, " = 0", "in a bivariate",
                    ifelse(type == "bvm", "sine von Mises", "wrapped Cauchy"))
    alt <- paste("H1:", dep, "!= 0")

  }

  # Construct an "htest" result
  result <- list(statistic = c("LRT statistic" = lr_stat), p.value = pvalue,
                 alternative = alt, method = method, df = df,
                 data.name = deparse(substitute(x)),
                 fit_mle = fit_mle, fit_null = fit_null)
  class(result) <- "htest"
  return(result)

}
