% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ridges-distr.R
\name{ridge_distr}
\alias{ridge_distr}
\alias{ridge_bvm}
\alias{ridge_bwc}
\alias{ridge_bwn}
\title{Connected component of the toroidal density ridge of a bivariate sine
von Mises, bivariate wrapped Cauchy, and bivariate wrapped normal}
\usage{
ridge_bvm(mu, kappa, eval_points, subint_1, subint_2)

ridge_bwc(mu, xi, eval_points, subint_1, subint_2)

ridge_bwn(mu, Sigma, kmax = 2, eval_points, subint_1, subint_2)
}
\arguments{
\item{mu}{circular means of the density, a vector of length \code{2}.}

\item{kappa}{vector of length \code{3} with the concentrations
\eqn{(\kappa_1, \kappa_2)} and the dependence parameter \eqn{\lambda}
of the density.}

\item{eval_points}{evaluation points for the ridge.}

\item{subint_1}{number of points for \eqn{\theta_1}.}

\item{subint_2}{number of points for \eqn{\theta_2} at each \eqn{\theta_1}.}

\item{xi}{a vector of length \code{3} with the marginal concentrations
\eqn{(\xi_1, \xi_2)}, and the dependence parameter \eqn{\rho}.}

\item{Sigma}{covariance matrix of size \code{c(2, 2)}.}

\item{kmax}{integer number up to truncate the wrapped normal series in
\code{-kmax:kmax}. Defaults to \code{2}.}
}
\value{
A matrix of size \code{c(subint_1, 2)} containing the points of the
 connected component of the ridge.
}
\description{
Computation of the connected component of the density ridge of
in a given set of points or, if not specified, in a regular grid on
\eqn{[-\pi, \pi)}.
}
\examples{
\donttest{
# Bivariate von Mises
mu <- c(0, 0)
kappa <- c(0.3, 0.5, 0.4)
nth <- 100
th <- seq(-pi, pi, l = nth)
x <- as.matrix(expand.grid(th, th))
d <- d_bvm(x = x, mu = mu, kappa = kappa)
image(th, th, matrix(d, nth, nth), col = viridisLite::viridis(20))
ridge <- ridge_bvm(mu = mu, kappa = kappa, subint_1 = 5e2,
                   subint_2 = 5e2)
points(ridge)

# Bivariate wrapped Cauchy
mu <- c(0, 0)
xi <- c(0.3, 0.6, 0.25)
nth <- 100
th <- seq(-pi, pi, l = nth)
x <- as.matrix(expand.grid(th, th))
d <- d_bwc(x = x, mu = mu, xi = xi)
image(th, th, matrix(d, nth, nth), col = viridisLite::viridis(20))
ridge <- ridge_bwc(mu = mu, xi = xi, subint_1 = 5e2, subint_2 = 5e2)
points(ridge)

# Bivariate wrapped normal
mu <- c(0, 0)
Sigma <- matrix(c(10, 3, 3, 5), nrow = 2)
nth <- 100
th <- seq(-pi, pi, l = nth)
x <- as.matrix(expand.grid(th, th))
d <- d_bwn(x = x, mu = mu, Sigma = Sigma)
image(th, th, matrix(d, nth, nth), col = viridisLite::viridis(20))
ridge <- ridge_bwn(mu = mu, Sigma = Sigma, subint_1 = 5e2,
                   subint_2 = 5e2)
points(ridge)}
}
\references{
Ozertem, U. and Erdogmus, D. (2011). Locally defined principal curves and
surfaces. \emph{Journal of Machine Learning Research}, 12(34):1249--1286.
\doi{10.6083/M4ZG6Q60}
}
