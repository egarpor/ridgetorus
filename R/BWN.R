
#' @title Density evaluation, sampling, and parameter estimation of the
#' bivariate wrapped normal distribution
#'
#' @description Computation of the density of the bivariate wrapped normal.
#'
#' @inheritParams bvm
#' @param Sigma covariance matrix of size \code{c(2, 2)}.
#' @param kmax integer number up to truncate the wrapped normal series in
#' \code{-kmax:kmax}. Defaults to \code{2}.
#' @return
#' \itemize{
#'   \item \code{d_bwn}: a vector of length \code{nx} with the density evaluated
#'   at \code{x}.
#'   \item \code{r_bwn}: a matrix of size \code{c(n, 2)} with the random sample.
#'   \item \code{fit_bwn_mle}: a list with the parameters
#'   \eqn{(\boldsymbol{\mu}, \boldsymbol{\Sigma})} and the object \code{opt}
#'   containing the optimization summary.
#' }
#' @examples
#' ## Density evaluation
#'
#' mu <- c(0, 0)
#' Sigma <- 3 * matrix(c(1.5, 0.75, 0.75, 1), nrow = 2, ncol = 2)
#' nth <- 50
#' th <- seq(-pi, pi, l = nth)
#' x <- as.matrix(expand.grid(th, th))
#' d <- d_bwn(x = x, mu = mu, Sigma = Sigma)
#' filled.contour(th, th, matrix(d, nth, nth), col = viridisLite::viridis(31),
#'                levels = seq(0, max(d), l = 30))
#' @export
#' @name bwn

#' @rdname bwn
#' @export
d_bwn <- function(x, mu, Sigma, kmax = 2) {

  # Center x
  x <- rbind(x)
  th1 <- sdetorus::toPiInt(x[, 1] - mu[1])
  th2 <- sdetorus::toPiInt(x[, 2] - mu[2])

  # Wrapping
  d <- numeric(length = nrow(x))
  k <- seq.int(-kmax, kmax)
  two_pi <- 2 * pi
  for (i in seq_along(k)) {
    for (j in seq_along(k)) {

      trans1 <- two_pi * k[i]
      trans2 <- two_pi * k[j]
      d <- d + mvtnorm::dmvnorm(x = cbind(th1 + trans1, th2 + trans2),
                                sigma = Sigma)

    }
  }
  return(d)

}


#' @description Simulation of pairs of samples from a bivariate wrapped normal.
#'
#' @inheritParams r_bvm
#' @rdname bwn
#' @export
r_bwn <- function(n, mu, Sigma) {

  sample <- mvtnorm::rmvnorm(n = n, mean = mu, sigma = Sigma)
  return(sdetorus::toPiInt(sample))

}


#' @description Maximum likelihood estimation of the parameters
#' \eqn{(\boldsymbol{\mu}, \boldsymbol{\Sigma})}.
#'
#' @inheritParams ridge_pca
#' @param lower,upper vector of length \code{5} with the bounds of the
#' parameters for the maximum likelihood optimizer.
#' @inheritParams fit_bvm_mle
#' @examples
#'
#' ## Sampling and estimation
#'
#' n <- 100
#' samp <- r_bwn(n = 100, mu = mu, Sigma = Sigma)
#' (param_mle <- fit_bwn_mle(samp)$par)
#' @rdname bwn
#' @export
fit_bwn_mle <- function(x, kmax = 2, lower = c(-pi, -pi, 1e-3, 1e-3, -1),
                        upper = c(pi, pi, 20, 20, 1), ...) {

  # Starting values
  ml1 <- circular::mle.wrappednormal(x = circular::circular(x[, 1]))
  ml2 <- circular::mle.wrappednormal(x = circular::circular(x[, 2]))
  start <- c(sdetorus::toPiInt(c(as.numeric(ml1$mu), as.numeric(ml2$mu))),
             ml1$sd^2, ml2$sd^2, 0)

  # Minus log-likelihood of the BWN
  minus_loglik_bwn <- function(param) {

    -sum(log(d_bwn(x = x, mu = param[1:2], kmax = kmax,
                   Sigma = matrix(c(param[3], param[5], param[5], param[4]),
                                  ncol = 2, nrow = 2))))

  }

  # Optimization
  opt <- sdetorus::mleOptimWrapper(minusLogLik = minus_loglik_bwn,
                                   start = start, lower = lower, upper = upper,
                                   optMethod = "Nelder-Mead", ...)

  return(list(mu = opt$par[1:2], Sigma = matrix(c(opt$par[3], opt$par[5],
                                                  opt$par[5], opt$par[4]),
                                                ncol = 2, nrow = 2)))

}
