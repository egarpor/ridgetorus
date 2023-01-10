
#' @title Density evaluation, sampling, and parameter estimation of the
#' bivariate wrapped Cauchy distribution
#'
#' @description Computation of the density of a bivariate wrapped Cauchy:
#' \deqn{f(\theta_1, \theta_2)=c(\xi_1,\xi_2,\rho)\{c_0(\xi_1,\xi_2,\rho)-
#' c_1(\xi_1,\xi_2,\rho) \cos (\theta_1-\mu_1)-
#' c_2(\xi_1,\xi_2,\rho)\cos (\theta_2-\mu_2)-\\
#' c_3(\xi_1,\xi_2,\rho) \cos (\theta_1-\mu_1) \cos (\theta_2-\mu_2)-
#' c_4(\xi_1,\xi_2,\rho) \sin (\theta_1-\mu_1) \sin (\theta_2-\mu_2)\}^{-1}.}
#'
#' @inheritParams bvm
#' @param xi a vector of length \code{3} with the marginal concentrations
#' \eqn{(\xi_1, \xi_2)}, and the dependence parameter \eqn{\rho}.
#' @return
#' \itemize{
#'   \item \code{d_bwc}: a vector of length \code{nx} with the density evaluated
#'   at \code{x}.
#'   \item \code{r_bwc}: a matrix of size \code{c(n, 2)} with the random sample.
#'   \item \code{fit_mme_bwc, fit_mle_bwc}: a list with the parameters
#'   \eqn{(\mu_1, \mu_2, \xi_1, \xi_2, \rho)} and the object \code{opt}
#'   containing the optimization summary.
#' }
#' @references
#' Kato, S. and Pewsey, A. (2015). A Möbius transformation-induced distribution
#' on the torus. \emph{Biometrika}, 102(2):359--370. \doi{10.1093/biomet/asv003}
#' @examples
#' ## Density evaluation
#'
#' mu <- c(0, 0)
#' xi <- c(0.3, 0.5, 0.4)
#' nth <- 50
#' th <- seq(-pi, pi, l = nth)
#' x <- as.matrix(expand.grid(th, th))
#' d <- d_bwc(x = x, mu = mu, xi = xi)
#' filled.contour(th, th, matrix(d, nth, nth), col = viridisLite::viridis(20),
#'                levels = seq(0, max(d), l = 20))
#' @name bwc


#' @rdname bwc
#' @export
d_bwc <- function(x, mu, xi) {

  # Extract parameters
  mu1 <- mu[1]
  mu2 <- mu[2]
  xi1 <- xi[1]
  xi2 <- xi[2]
  rho <- xi[3]

  # Constants
  rho2 <- rho^2
  xi12 <- xi1^2
  xi22 <- xi2^2
  C <- (1 - rho2) * (1 - xi12) * (1 - xi22) / (4 * pi^2)
  c0  <-  (1 + rho2)  *  (1 + xi12)  *  (1 + xi22)  -
    8  *  abs(rho) * xi1 * xi2
  c1 <- 2 * (1 + rho2) * xi1 * (1 + xi22) - 4 * abs(rho) * (1 + xi12) * xi2
  c2 <- 2 * (1 + rho2) * (1 + xi12) * xi2 - 4 * abs(rho) * xi1 * (1 + xi22)
  c3 <- -4 * (1 + rho2) * xi1 * xi2 + 2 * abs(rho) * (1 + xi12) * (1 + xi22)
  c4 <- 2 * rho * (1 - xi12) * (1 - xi22)

  # Density
  x <- rbind(x)
  cos1 <- cos(x[, 1] - mu1)
  cos2 <- cos(x[, 2] - mu2)
  return(C / (c0 - c1 * cos1 - c2 * cos2 - c3 * cos1 * cos2 -
                c4 * sin(x[, 1] - mu1) * sin(x[, 2] - mu2)))

}


#' @description Simulation of samples from a bivariate wrapped Cauchy.
#'
#' @inheritParams r_bvm
#' @author The original code for \code{r_bwc} was supplied by
#' Arthur Pewsey.
#' @rdname bwc
#' @export
r_bwc <- function(n, mu, xi) {

  # Uses Theorem 1 to simulate a random bivariate complex vector (Z1, Z2)

  # Simulate Z1 with marginal density and equation (8)
  xi1 <- xi[1]
  xi2 <- xi[2]
  rho <- xi[3]
  phi <- abs(rho)
  U1 <- runif(n, 0, 1)
  C1 <- exp(1i * 2 * pi * U1)
  Z1 <- (C1 + xi1) / (xi1 * C1 + 1)

  # Obtain the first circular variable as its argument
  theta1 <- sdetorus::toPiInt(Arg(Z1) + mu[1])

  # Simulate Z2 using the conditional marginal and equation (8)
  C1inv <- 1 / C1
  tmp <- matrix(runif(2 * n, 0, 1), ncol = 2, byrow = TRUE)
  U2 <- tmp[, 1]
  C2 <- exp(1i * 2 * pi * U2)
  if (rho < 0) {

    alpha2 <- (phi * xi2 * C1 + 1) / (phi * xi2 * C1inv + 1)
    beta2 <- (xi2 + phi * C1inv) / (phi * xi2 * C1 + 1)

  } else if (rho >= 0) {

    alpha2 <- (phi * xi2 * C1inv + 1) / (phi * xi2 * C1 + 1)
    beta2 <- (xi2 + phi * C1) / (phi * xi2 * C1inv + 1)

  }
  Z2 <- alpha2 * (C2 + beta2) / (Conj(beta2) * C2 + 1)

  # Obtain the second circular variable as its argument
  theta2 <- sdetorus::toPiInt(Arg(Z2) + mu[2])

  # Return the pair of circular variables
  return(cbind(theta1, theta2))

}


#' @description Maximum likelihood and method of moments estimation of the
#' parameters \eqn{(\mu_1, \mu_2, \xi_1, \xi_2, \rho)}.
#'
#' @inheritParams ridge_pca
#' @param start a vector of length \code{5} with the initial values for the
#' maximum likelihood optimizer. If \code{NULL} (default), the method of
#' moments estimates are employed.
#' @param lower,upper vectors of length \code{5} with the bounds for the
#' likelihood optimizer. Default to \code{c(-pi, -pi, 0, 0, -1 + 1e-3)} and
#' \code{c(pi, pi, 1 - 1e-3, 1 - 1e-3, 1 - 1e-3)}.
#' @inheritParams bvm
#' @examples
#'
#' ## Sampling and estimation
#'
#' n <- 100
#' samp <- r_bwc(n = n, mu = mu, xi = xi)
#' (param_mm <- fit_bwc_mm(samp)$par)
#' (param_mle <- fit_bwc_mle(samp)$par)
#' @rdname bwc
#' @export
fit_bwc_mm <- function(x, hom = FALSE, indep = FALSE) {

  # Circular means
  n <- nrow(x)
  theta1 <- x[, 1]
  theta2 <- x[, 2]
  R1_cos <- sum(cos(theta1))
  R1_sin <- sum(sin(theta1))
  R2_cos <- sum(cos(theta2))
  R2_sin <- sum(sin(theta2))
  mu1_hat <- atan2(x = R1_cos, y = R1_sin)
  mu2_hat <- atan2(x = R2_cos, y = R2_sin)

  # Concentration parameters
  if (hom) {

    xi1_hat <- sqrt((R1_cos + R2_cos)^2 + (R1_sin + R2_sin)^2) / (2 * n)
    xi2_hat <- xi1_hat

  } else {

    xi1_hat <- sqrt(R1_cos^2 + R1_sin^2) / n
    xi2_hat <- sqrt(R2_cos^2 + R2_sin^2) / n

  }

  # Dependence parameter
  if (indep) {

    rho_hat <- 0

  } else {

    phi1 <- 2 * atan((1 + xi1_hat) / (1 + xi1_hat) *
                       tan((theta1 - mu1_hat) / 2))
    phi2 <- 2 * atan((1 + xi2_hat) / (1 + xi2_hat) *
                       tan((theta2 - mu2_hat) / 2))
    rho_hat <- (sqrt(sum(cos(phi1 - phi2))^2 + sum(sin(phi1 - phi2))^2) -
                  sqrt(sum(cos(phi1 + phi2))^2 + sum(sin(phi1 + phi2))^2)) / n

  }

  # Result
  return(list(mu1 = mu1_hat, mu2 = mu2_hat, xi1 = xi1_hat, xi2 = xi2_hat,
              rho = rho_hat))

}


#' @rdname bwc
#' @export
fit_bwc_mle <- function(x, start = NULL, lower = c(-pi, -pi, 0, 0, -1 + 1e-3),
                        upper = c(pi, pi, 1 - 1e-3, 1 - 1e-3, 1 - 1e-3),
                        hom = FALSE, indep = FALSE, ...) {

  # Minus log-likelihood of the BWC
  minus_loglik_bwc <- function(param) {

    # Set parameters. Homogeneous concentrations? Independence?
    mu <- param[1:2]
    xi1 <- param[3]
    xi2 <- ifelse(hom, xi1, param[4])
    rho <- ifelse(indep, 0, param[length(param)])

    # Minus log-likelihood
    return(-sum(log(d_bwc(x = x, mu = mu, xi = c(xi1, xi2, rho)))))

  }

  # Fits if independence
  if (indep) {

    # Marginal fits
    mle_1 <- circular:::MlewrappedcauchyRad(x = x[, 1], mu = NULL, rho = NULL,
                                            tol = 1e-10, max.iter = 200)
    mle_2 <- circular:::MlewrappedcauchyRad(x = x[, 2], mu = NULL, rho = NULL,
                                            tol = 1e-10, max.iter = 200)
    param <- c(sdetorus::toPiInt(c(mle_1[1], mle_2[1])),
               mle_1[2], mle_2[2], 0)

    # Update parameters if homogeneity
    if (hom) {

      # Reduce start/lower/upper if homogeneity
      start <- c(param[1:2], mean(param[3:4]))
      lower <- c(lower[1:2], mean(lower[3:4]))
      upper <- c(upper[1:2], mean(upper[3:4]))

      # Optimization
      opt <- sdetorus::mleOptimWrapper(minusLogLik = minus_loglik_bwc,
                                       optMethod = "L-BFGS-B", start = start,
                                       selectSolution = "lowest", lower = lower,
                                       upper = upper, checkCircular = TRUE, ...)
      param <- opt$par

    }

    # Fitted parameters
    mu1 <- param[1]
    mu2 <- param[2]
    xi1 <- param[3]
    xi2 <- ifelse(hom, xi1, param[4])
    rho <- 0

    # Optimization function
    opt <- list(par = param, value = minus_loglik_bwc(param = param))

  } else {

    # If not provided, obtain starting values by moment estimators
    x <- rbind(x)
    if (is.null(start)) {

      start <- unlist(fit_bwc_mm(x = x, hom = hom, indep = indep))

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
    opt <- sdetorus::mleOptimWrapper(minusLogLik = minus_loglik_bwc,
                                     optMethod = "L-BFGS-B", start = start,
                                     selectSolution = "lowest", lower = lower,
                                     upper = upper, checkCircular = TRUE, ...)

    # Fitted parameters
    mu1 <- opt$par[1]
    mu2 <- opt$par[2]
    xi1 <- opt$par[3]
    xi2 <- ifelse(hom, xi1, opt$par[4])
    rho <- opt$par[length(opt$par)]

  }

  # Result
  return(list(mu1 = mu1, mu2 = mu2, xi1 = xi1, xi2 = xi2, rho = rho, opt = opt))

}
