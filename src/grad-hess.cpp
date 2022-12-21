
// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::plugins("cpp11")]]
#include <RcppArmadillo.h>
#include <Rmath.h>
#include <R.h>

using namespace Rcpp;
using namespace std;

//' @title Unnormalized gradient and Hessian of a multivariate sine von Mises
//'
//' @description Computation of the gradient and Hessian of an arbitrary
//' multivariate sine von Mises density, without normalizing constants and
//' density factor.
//'
//' @param theta a matrix of size \code{c(nx, d)} with angles on
//' \eqn{[-\pi, \pi)}.
//' @param kappa vector with the \eqn{d} concentration parameters
//' \eqn{\boldsymbol{\kappa} = (\kappa_1, \ldots, \kappa_d)'}.
//' @param Lambda dependence matrix \eqn{\boldsymbol{\Lambda}}.
//' @return A list:
//' \item{grad}{unnormalized gradient, a matrix of size \code{c(nx, d)}.}
//' \item{hess}{unnormalized Hessian, an array of size \code{c(nx, d, d)}.}
//' @references
//' Mardia, K. V., Hughes, G., Taylor, C. C., and Singh, H. (2008).
//' A multivariate von Mises with applications to bioinformatics.
//' \emph{Canadian Journal of Statistics}, 36(1):99--109.
//' \doi{10.1002/cjs.5550360110}
//' @noRd
// [[Rcpp::export]]
Rcpp::List grad_hess_mvm(arma::mat theta, arma::vec kappa, arma::mat Lambda) {

  // Dimensions
  arma::uword nx = theta.n_rows;
  arma::uword d = theta.n_cols;

  // Needed objects
  arma::mat delta = arma::eye(d, d);
  arma::mat c = cos(theta);
  arma::mat s = sin(theta);

  // Without normalization constant and f_mvm common factor
  arma::mat grad(nx, d);
  arma::cube hess(nx, d, d);
  for(arma::uword i = 0; i < d; i++) {

    // Gradient
    arma::vec sum = (Lambda.row(i) * s.t()).t();
    grad.col(i) = -kappa[i] * s.col(i) + c.col(i) % sum;

    // Hessian
    for(arma::uword j = 0; j < d; j++) {

      hess(arma::span::all, arma::span(i), arma::span(j)) =
        (-kappa[i] * c.col(i) + s.col(i) % sum) * delta(i, j) +
        c.col(i) % c.col(j) * Lambda(i, j);

    }
  }

  return Rcpp::List::create(Rcpp::Named("grad") = grad,
                            Rcpp::Named("hess") = hess);

}


//' @title Unnormalized gradient and Hessian of a multivariate wrapped normal
//'
//' @description Computation of the gradient and Hessian of an arbitrary
//' multivariate wrapped normal density, without normalizing constants and
//' density factor.
//'
//' @inheritParams grad_hess_mvm
//' @param mu vector of length \code{d} with the mean of the normal
//' distribution.
//' @param Sigma matrix of size \code{c(d, d)} with the covariance matrix of
//' the normal distribution.
//' @param k integer values for the wrapped normal truncation.
//' @return A list:
//' \item{grad}{unnormalized gradient, a matrix of size \code{c(nx, d)}.}
//' \item{hess}{unnormalized Hessian, an array of size \code{c(nx, d, d)}.}
//' @noRd
// [[Rcpp::export]]
Rcpp::List grad_hess_mwn(arma::mat theta, arma::vec mu, arma::mat Sigma,
                         arma::mat k) {

  // Center evaluation points with mu
  arma::mat inv_sigma = inv(Sigma);
  arma::uword nx = theta.n_rows;
  arma::uword d = mu.n_elem;
  arma::uword nrow = k.n_rows;
  arma::mat theta_mu(nx, d);
  for (arma::uword i = 0; i < d; i++) {

    theta_mu.col(i) = theta.col(i) - mu[i];

  }

  // Results
  arma::mat grad(nx, d);
  arma::cube hess(nx, d, d);

  // Loop on observations
  for(arma::uword i = 0; i < nx; i++) {

    // Approximate the result as the sum for k = {-1, 0, 1}
    arma::vec grad_aux = arma::zeros(d);
    arma::mat hess_aux = arma::zeros(d, d);
    for(arma::uword m = 0; m < nrow; m++) {

      // Factor 1 / ((2 * pi) |Sigma|) not taken into account
      arma::mat theta_f = 2 * M_PI * k.row(m) + theta_mu.row(i);
      double dens = as_scalar(exp(-0.5 * theta_f * inv_sigma * theta_f.t()));
      grad_aux -= dens * inv_sigma * theta_f.t();
      hess_aux += dens * inv_sigma * (theta_f.t() * theta_f - Sigma) *
        inv_sigma;

    }
    grad.row(i) = grad_aux.t();
    hess(arma::span(i), arma::span::all, arma::span::all) = hess_aux;

  }

  return Rcpp::List::create(Rcpp::Named("grad") = grad,
                            Rcpp::Named("hess") = hess);

}


//' @title Unnormalized gradient and Hessian of a bivariate wrapped Cauchy
//'
//' @description Gradient and Hessian of a bivariate wrapped Cauchy density,
//' without normalizing constants and density factor.
//'
//' @param theta2 evaluation points (vector).
//' @param theta1 evaluation point (scalar).
//' @inheritParams bwc
//' @return A list:
//' \item{D1,D2}{vectors of size \code{nx} with the unnormalized
//' first derivative.}
//' \item{u,w}{vector of size \code{nx} with the unnormalized
//' second derivative.}
//' \item{v}{vector of size \code{nx} with the unnormalized mixed derivative.}
//' @noRd
// [[Rcpp::export]]
Rcpp::List grad_hess_bwc(arma::vec theta2, double theta1, arma::vec xi) {

  // Compute constants
  double rho2 = pow(xi(2), 2);
  double xi12 = pow(xi(0), 2);
  double xi22 = pow(xi(1), 2);
  double absrho = abs(xi(2));
  double c0 = (1 + rho2) * (1 + xi12) * (1 + xi22) -
   8 * absrho * xi(0) * xi(1);
  double c1 = 2 * (1 + rho2) * xi(0) * (1 + xi22) -
   4 * absrho * (1 + xi12) * xi(1);
  double c2 = 2 * (1 + rho2) * (1 + xi12) * xi(1) -
   4 * absrho * xi(0) * (1 + xi22);
  double c3 = -4 * (1 + rho2) * xi(0) * xi(1) +
   2 * absrho * (1 + xi12) * (1 + xi22);
  double c4 = 2 * xi(2) * (1 - xi12) * (1 - xi22);

  // Inverse density
  double sin1 = -sin(theta1);
  double cos1 = cos(theta1);
  arma::vec sin2 = -sin(theta2);
  arma::vec cos2 = cos(theta2);
  arma::vec sin12 = sin1 * sin2;
  arma::vec cos12 = cos1 * cos2;
  arma::vec sin1cos2 = sin1 * cos2;
  arma::vec sin2cos1 = sin2 * cos1;
  arma::vec d = c0 - c1 * cos1 - c2 * cos2 - c3 * cos12 - c4 * sin12;

  // Gradients
  arma::vec D1 = (c1 * sin1 + c3 * sin1cos2 - c4 * sin2cos1);
  arma::vec D2 = (c2 * sin2 + c3 * sin2cos1 - c4 * sin1cos2);

  // Second derivatives
  arma::vec u = 2 * (c1 * sin1 + c3 * sin1cos2 - c4 * sin2cos1) %
   (c1 * sin1 + c3 * sin1cos2 - c4 * sin2cos1) / arma::pow(d, 3) +
   (-c1 * cos1 - c3 * cos12 - c4 * sin12) / arma::square(d);
  arma::vec w = 2 * (c2 * sin2 + c3 * sin2cos1 - c4 * sin1cos2) %
   (c2 * sin2 + c3 * sin2cos1 - c4 * sin1cos2) / arma::pow(d, 3) +
   (-c2 * cos2 - c3 * cos12 - c4 * sin12) / arma::square(d);
  arma::vec v = (c3 * sin12 + c4 * cos12) / arma::square(d) +
   2 * (c1 * sin1 + c3 * sin1cos2 - c4 * sin2cos1) %
   (c2 * sin2 + c3 * sin2cos1 - c4 * sin1cos2) / arma::pow(d, 3);

  return Rcpp::List::create(Rcpp::Named("D1") = D1, Rcpp::Named("D2") = D2,
                            Rcpp::Named("u") = u, Rcpp::Named("v") = v,
                            Rcpp::Named("w") = w);

}
