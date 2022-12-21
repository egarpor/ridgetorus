
#' @title Density evaluation, sampling, and parameter estimation of the
#' bivariate sine von Mises distribution
#'
#' @description Computation of the density and normalizing constant
#' \eqn{T(\kappa_1, \kappa_2, \lambda)} of the bivariate sine von Mises
#' \deqn{f(\theta_1, \theta_2)= T(\kappa_1, \kappa_2, \lambda)
#' \exp\{\kappa_1 \cos(\theta_1-\mu_1) +
#' \kappa_2 \cos(\theta_2-\mu_2) +
#' \lambda \sin(\theta_1-\mu_1) \sin(\theta_2-\mu_2)\}.}
#'
#' @param x matrix of size \code{c(nx, 2)} with the angles on which the density
#' is evaluated.
#' @param mu circular means of the density, a vector of length \code{2}.
#' @param kappa vector of length \code{3} with the concentrations
#' \eqn{(\kappa_1, \kappa_2)} and the dependence parameter \eqn{\lambda}
#' of the density.
#' @param log_const logarithm of the normalizing constant. Computed internally
#' if \code{NULL} (default).
#' @param M truncation of the series expansion for computing the normalizing
#' constant. Defaults to \code{25}.
#' @param MC Monte Carlo replicates for computing the normalizing
#' constant when there is no series expansion. Defaults to \code{1e4}.
#' @return
#' \itemize{
#'   \item \code{d_bvm}: a vector of length \code{nx} with the density evaluated
#'   at \code{x}.
#'   \item \code{const_bvm}: the value of the normalizing constant
#'   \eqn{T(\kappa_1, \kappa_2, \lambda)}.
#'   \item \code{r_bvm}: a matrix of size \code{c(n, 2)} with the random sample.
#'   \item \code{fit_mme_bvm, fit_mle_bvm}: a list with the estimated parameters
#'   \eqn{(\mu_1, \mu_2, \kappa_1, \kappa_2, \lambda)} and the object \code{opt}
#'   containing the optimization summary.
#' }
#' @references
#' Mardia, K. V., Hughes, G., Taylor, C. C., and Singh, H. (2008).
#' A multivariate von Mises with applications to bioinformatics.
#' \emph{Canadian Journal of Statistics}, 36(1):99--109.
#' \doi{10.1002/cjs.5550360110}
#'
#' Singh, H., Hnizdo, V., and Demchuk, E. (2002). Probabilistic model for two
#' dependent circular variables. \emph{Biometrika}, 89(3):719--723.
#' \doi{10.1093/biomet/89.3.719}
#' @examples
#' ## Density evaluation
#'
#' mu <- c(0, 0)
#' kappa <- 3:1
#' nth <- 50
#' th <- seq(-pi, pi, l = nth)
#' x <- as.matrix(expand.grid(th, th))
#' const <- const_bvm(kappa = kappa)
#' d <- d_bvm(x = x, mu = mu, kappa = kappa, log_const = log(const))
#' filled.contour(th, th, matrix(d, nth, nth), col = viridisLite::viridis(31),
#'                levels = seq(0, max(d), l = 30))
#' @name bvm


#' @rdname bvm
#' @export
d_bvm <- function(x, mu, kappa, log_const = NULL) {

  # As matrix
  x <- rbind(x)

  # Normalizing constant
  if (is.null(log_const)) {

    log_const <- log(const_bvm(kappa = kappa))

  }

  # Differences
  th_mu_1 <- x[, 1] - mu[1]
  th_mu_2 <- x[, 2] - mu[2]

  # Density
  exp(log_const + kappa[1] * cos(th_mu_1) + kappa[2] * cos(th_mu_2) +
        kappa[3] * sin(th_mu_1) * sin(th_mu_2))

}


#' @rdname bvm
#' @export
const_bvm <- function(kappa, M = 25, MC = 1e4) {

  # Cases: independence; safe series expansion; Monte Carlo.
  m <- 0:M
  if (abs(kappa[3]) < 1e-10) {

    const <- 4 * pi^2 * exp(sum(
      kappa[1:2] + log(besselI(x = kappa[1:2], nu = 0, expon.scaled = TRUE))))

  } else if (abs(kappa[3]) <= 30) {

    kappa[kappa == 0] <- 1e-10
    const <- 4 * pi^2 * sum(exp(lchoose(n = 2 * m, k = m) +
      m * (2 * (log(abs(kappa[3])) - log(2)) - log(kappa[1]) - log(kappa[2])) +
      kappa[1] + log(besselI(x = kappa[1], nu = m, expon.scaled = TRUE)) +
      kappa[2] + log(besselI(x = kappa[2], nu = m, expon.scaled = TRUE))))

  } else {

    message(paste("Unreliable constant series for kappa1 = 0 or kappa2 = 0 and",
                  "abs(lambda) <= 30. Using Monte Carlo integration."))
    const <- sdetorus::mcTorusIntegrate(f = function(x) {
      d_bvm(x = x, mu = c(0, 0), kappa = kappa, log_const = 0)},
      p = 2, M = MC, fVect = TRUE)

  }

  # Return the actual normalizing constant
  return(1 / const)

}


#' @description Simulation of samples from a bivariate sine von Mises.
#'
#' @param n sample size.
#' @rdname bvm
#' @export
r_bvm <- function(n, mu, kappa) {

  sample <- BAMBI::rvmsin(n = n, kappa1 = kappa[1], kappa2 = kappa[2],
                          kappa3 = kappa[3], mu1 = mu[1], mu2 = mu[2])
  return(sdetorus::toPiInt(sample))

}


#' @description Maximum likelihood and method of moments estimation of the
#' parameters \eqn{(\mu_1, \mu_2, \kappa_1, \kappa_2, \lambda)}.
#'
#' @param start a vector of length \code{5} with the initial values for the
#' maximum likelihood optimizer. The first two entries are disregarded in
#' \code{fit_bvm_mm}. If \code{NULL} (default), the starting values are taken
#' from the estimation of marginal von Mises in \code{fit_bvm_mm}. In
#' \code{fit_bvm_mle}, the method of moments estimates are employed.
#' @param lower,upper vectors of length \code{5} with the bounds for the
#' likelihood optimizer. Default to \code{c(-pi, -pi, 0, 0, -30)} and
#' \code{c(pi, pi, 30, 30, 30)}.
#' @param ... further parameters passed to
#' \code{\link[=sdetorus]{mleOptimWrapper}}.
#' @param hom assume a homogeneous distribution with equal marginal
#' concentrations? Defaults to \code{FALSE}.
#' @param indep set the dependence parameter to zero? Defaults to \code{FALSE}.
#' @examples
#'
#' ## Sampling and estimation
#'
#' n <- 100
#' samp <- r_bvm(n = n, mu = mu, kappa = kappa)
#' (param_mm <- fit_bvm_mm(samp)$par)
#' (param_mle <- fit_bvm_mle(samp)$par)
#' @rdname bvm
#' @export
fit_bvm_mm <- function(x, lower = c(0, 0, -30), upper = c(30, 30, 30),
                       start = NULL, M = 25, hom = FALSE, indep = FALSE, ...) {

  # Sensible start values
  if (is.null(start)) {

    start <- c(DirStats::kappa_ml(DirStats::to_cir(x[, 1])),
               DirStats::kappa_ml(DirStats::to_cir(x[, 2])), 0)

  } else {

    l_start <- length(start)
    if (l_start == 5) {

      start <- start[-c(1:2)]

    } else if (l_start != 3) {

      stop(paste("start must have length 3 (concentrations) or",
                 "5 (if means are included, ignored)."))

    }

  }

  # Alter start values if homogeneity or independence
  if (hom) {

    start <- c(mean(start[1:2]), start[3])
    lower <- c(mean(lower[1:2]), lower[3])
    upper <- c(mean(upper[1:2]), upper[3])

  }
  if (indep) {

    start <- start[-length(start)]
    lower <- lower[-length(lower)]
    upper <- upper[-length(upper)]

  }

  # Circular means
  R1_cos <- sum(cos(x[, 1]))
  R1_sin <- sum(sin(x[, 1]))
  R2_cos <- sum(cos(x[, 2]))
  R2_sin <- sum(sin(x[, 2]))
  mu1 <- atan2(x = R1_cos, y = R1_sin)
  mu2 <- atan2(x = R2_cos, y = R2_sin)

  # Center data
  th1 <- x[, 1] - mu1
  th2 <- x[, 2] - mu2

  # We apply equations (10)-(12) in Mardia et al. (2008), draft version from
  # https://citeseerx.ist.psu.edu/pdf/63b0c8c5c98acab0d8866af8518470306363bc64
  # Note lambda_12 is 0.5 * lambda, there is an inconsistency on the moment
  # equations and equation (1) and the following equation to (12).

  # LHS
  cos1 <- mean(cos(th1))
  cos2 <- mean(cos(th2))
  sin12 <- mean(sin(th1) * sin(th2))

  # Moment equations for concentrations
  squared_moments_bvm <- function(param) {

    # Independence? Homogeneous concentrations?
    kappa1 <- param[1]
    kappa2 <- ifelse(hom, param[1], param[2])
    lambda <- ifelse(indep, 0, param[length(param)])

    # Common part in series
    m <- 0:(M + 1)
    log_bes_1 <- kappa1 + log(besselI(x = kappa1, nu = m, expon.scaled = TRUE))
    log_bes_2 <- kappa2 + log(besselI(x = kappa2, nu = m, expon.scaled = TRUE))
    c <- c(0, lchoose(n = 2 * m[2:(M + 1)], k = m[2:(M + 1)]) +
             m[2:(M + 1)] * (2 * log(abs(0.5 * lambda)) -
                               (log(kappa1) + log(kappa2))))
    den <- sum(exp(c[1:M] + log_bes_1[1:M] + log_bes_2[1:M]), na.rm = TRUE)

    # RHS
    rhs1 <- sum(exp(c[1:M] + log_bes_1[2:(M + 1)] + log_bes_2[1:M]),
                na.rm = TRUE) / den
    rhs2 <- sum(exp(c[1:M] + log_bes_1[1:M] + log_bes_2[2:(M + 1)]),
                na.rm = TRUE) / den
    rhs3 <- sum(exp(log(m[2:(M + 1)]) + c[2:(M + 1)] +
                      log_bes_1[2:(M + 1)] + log_bes_2[2:(M + 1)]),
                na.rm = TRUE) / (0.5 * lambda * den)

    # Return sum(lhs - rhs)^2
    return((cos1 - rhs1)^2 + (cos2 - rhs2)^2 + (sin12 - rhs3)^2)

  }

  # Optimization
  opt <- sdetorus::mleOptimWrapper(minusLogLik = squared_moments_bvm,
                                   optMethod = "L-BFGS-B", start = start,
                                   selectSolution = "lowest", lower = lower,
                                   upper = upper, checkCircular = TRUE, ...)

  # Enforce tweaked parameters
  kappa1 <- opt$par[1]
  kappa2 <- ifelse(hom, kappa1, opt$par[2])
  lambda <- ifelse(indep, 0, opt$par[length(opt$par)])

  # Result
  return(list(mu1 = mu1, mu2 = mu2, kappa1 = kappa1, kappa2 = kappa2,
              lambda = lambda))

}


#' @rdname bvm
#' @export
fit_bvm_mle <- function(x, start = NULL, M = 25, lower = c(-pi, -pi, 0, 0, -30),
                        upper = c(pi, pi, 30, 30, 30), hom = FALSE,
                        indep = FALSE, ...) {

  # Minus log-likelihood of the BvM
  minus_loglik_bvm <- function(param) {

    # Set parameters. Homogeneous concentrations? Independence?
    mu <- param[1:2]
    kappa1 <- param[3]
    kappa2 <- ifelse(hom, kappa1, param[4])
    lambda <- ifelse(indep, 0, param[length(param)])

    # Minus log-likelihood
    log_const <- log(const_bvm(kappa = c(kappa1, kappa2, lambda), M = M))
    return(-sum(log(d_bvm(x = x, mu = mu, kappa = c(kappa1, kappa2, lambda),
                          log_const = log_const))))

  }

  # Fits if independence
  if (indep) {

    lambda <- 0
    x1 <- DirStats::to_cir(th = x[, 1])
    x2 <- DirStats::to_cir(th = x[, 2])
    mu1 <- sdetorus::toPiInt(DirStats::to_rad(DirStats::mu_ml(x1)))
    mu2 <- sdetorus::toPiInt(DirStats::to_rad(DirStats::mu_ml(x2)))

    if (hom) {

      R1 <- DirStats::norm2(colMeans(x1))
      R2 <- DirStats::norm2(colMeans(x2))
      R <- (R1 + R2) / 2
      A_2 <- function(kappa) besselI(x = kappa, nu = 1, expon.scaled = TRUE) /
          besselI(x = kappa, nu = 0, expon.scaled = TRUE) - R
      kappa1 <- uniroot(A_2, lower = 1e-04, upper = 100, tol = 1e-04)$root
      kappa2 <- kappa1
      param <- c(mu1, mu2, kappa1, lambda)

    } else {

      kappa1 <- DirStats::kappa_ml(x1)
      kappa2 <- DirStats::kappa_ml(x2)
      param <- c(mu1, mu2, kappa1, kappa2, lambda)

    }

    # Optimization function
    opt <- list(par = param, value = minus_loglik_bvm(param = param))

  } else {

    # If not provided, obtain starting values by moment estimators
    x <- rbind(x)
    if (is.null(start)) {

      start <- unlist(fit_bvm_mm(x = x, hom = hom, indep = FALSE))

    } else {

      if (length(start) != 5) {

        stop("start must have length 5.")

      }

    }

    # Reduce start/lower/upper if homogeneity
    if (hom) {

      start <- c(start[1:2], mean(start[3:4]), start[5])
      lower <- c(lower[1:2], mean(lower[3:4]), lower[5])
      upper <- c(upper[1:2], mean(upper[3:4]), upper[5])

    }

    # Optimization
    opt <- sdetorus::mleOptimWrapper(minusLogLik = minus_loglik_bvm,
                                     optMethod = "L-BFGS-B", start = start,
                                     selectSolution = "lowest", lower = lower,
                                     upper = upper, checkCircular = TRUE, ...)

    # Fitted parameters
    mu1 <- opt$par[1]
    mu2 <- opt$par[2]
    kappa1 <- opt$par[3]
    kappa2 <- ifelse(hom, kappa1, opt$par[4])
    lambda <- opt$par[length(opt$par)]

  }

  # Result
  return(list(mu1 = mu1, mu2 = mu2, kappa1 = kappa1, kappa2 = kappa2,
              lambda = lambda, opt = opt))

}
