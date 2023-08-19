
// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <Rmath.h>
#include <R.h>

using namespace Rcpp;
using namespace std;

// Declaration for functions
Rcpp::List grad_hess_mvm(arma::mat theta, arma::vec kappa, arma::mat Lambda);
Rcpp::List grad_hess_bwc(arma::vec theta2, double theta1, arma::vec xi);
Rcpp::List grad_hess_mwn(arma::mat theta, arma::vec mu, arma::mat Sigma,
                         arma::mat k);


//' @title Implicit equation of a toroidal density ridge
//'
//' @description One of the conditions for the density ridge of a given density
//' \eqn{f} is that the modulus of its projected gradient is zero:
//' \eqn{\|\mathrm{D}_{p-1}f(\mathbf{x})\|=0}. This function computes the
//' LHS of that implicit equation for the case of a given density.
//'
//' @param theta1,theta2 evaluation points.
//' @param density chosen model for the density: \code{"bvm"}, \code{"bwc"}, or
//' \code{"bwn"}.
//' @inheritParams d_bwc
//' @param kappa,Lambda vector of concentrations and dependence matrix of the
//' multivariate von Mises distribution.
//' @param mu,Sigma vector of means and covariance matrix of the multivariate
//' normal distribution.
//' @inheritParams grad_hess_mwn
//' @return The value of the LHS of the implicit equation.
//' @examples
//' n <- 200
//' x <- seq(-pi, pi, l = n)
//' mu <- c(0, 0)
//' kappa <- c(0.3, 0.4, 0.5)
//' val <- sapply(x, function(th1) ridgetorus:::implicit_equation(
//'     theta2 = x, theta1 = th1, density = "bvm", kappa = kappa[1:2],
//'     Lambda = matrix(c(0, kappa[3], kappa[3], 0), nrow = 2, ncol = 2)))
//' val <- matrix(val, nrow = n, ncol = n)
//' old_par <- par(no.readonly = TRUE)
//' par(mfrow = c(1, 2))
//' image(x, x, -log(abs(val)), axes = FALSE, col = viridisLite::viridis(20))
//' sdetorus::torusAxis()
//' sdetorus::plotSurface2D(x, x, f = function(x) d_bvm(x = x, mu = mu,
//'                                                     kappa = kappa),
//'                         axes = FALSE)
//' sdetorus::torusAxis()
//' par(old_par)
//' @noRd
// [[Rcpp::export]]
arma::vec implicit_equation(arma::vec theta2, double theta1, String density,
                            arma::vec kappa = 0,
                            Rcpp::Nullable<NumericMatrix> Lambda = R_NilValue,
                            arma::vec xi = 0, arma::vec mu = 0,
                            Rcpp::Nullable<NumericMatrix> Sigma = R_NilValue,
                            Rcpp::Nullable<NumericMatrix> k = R_NilValue) {

  // Stop if theta2 is empty
  if (theta2.n_elem == 0) {

    Rcpp::stop("theta2 is empty");

  }

  // Needed declarations
  Rcpp::List gradhess;
  arma::vec D1 = arma::zeros(theta2.n_elem);
  arma::vec D2 = arma::zeros(theta2.n_elem);
  arma::vec u = arma::zeros(theta2.n_elem);
  arma::vec v = arma::zeros(theta2.n_elem);
  arma::vec w = arma::zeros(theta2.n_elem);

  if (density == "bvm" || density == "bwn") { // Valid for multivariate

    arma::mat theta = arma::zeros(theta2.n_elem, 2);
    theta.col(0).fill(theta1);
    theta.col(1) = theta2;
    if (density == "bvm") {

      gradhess = grad_hess_mvm(theta, kappa,
                               Rcpp::as<arma::mat>(Rcpp::wrap(Lambda)));

    } else {

      gradhess = grad_hess_mwn(theta, mu,
                               Rcpp::as<arma::mat>(Rcpp::wrap(Sigma)),
                               Rcpp::as<arma::mat>(Rcpp::wrap(k)));

    }
    arma::mat grad = gradhess["grad"];
    arma::cube hess = gradhess["hess"];
    D1 = grad.col(0);
    D2 = grad.col(1);
    u = hess(arma::span::all, arma::span(0), arma::span(0));
    v = hess(arma::span::all, arma::span(0), arma::span(1));
    w = hess(arma::span::all, arma::span(1), arma::span(1));

  } else if (density == "bwc") { // Different because it is only bivariate

    gradhess = grad_hess_bwc(theta2, theta1, xi);
    D1 = Rcpp::as<arma::vec>(gradhess["D1"]);
    D2 = Rcpp::as<arma::vec>(gradhess["D2"]);
    u = Rcpp::as<arma::vec>(gradhess["u"]);
    v = Rcpp::as<arma::vec>(gradhess["v"]);
    w = Rcpp::as<arma::vec>(gradhess["w"]);

  } else {

    Rcpp::stop("Unknown model");

  }

  // Eigenvectors
  arma::vec G1 = 2 * (u - w + v - sqrt(arma::square(w - u) +
    4 * arma::square(v)));
  arma::vec G2 = w - u + 4 * v - sqrt(arma::square(w - u) +
    4 * arma::square(v));

  // Normalize
  arma::vec G_norm = sqrt(arma::square(G1) + arma::square(G2));
  G1 /= G_norm;
  G2 /= G_norm;

  // Implicit ridge expression
  return(D1 % G1 + D2 % G2);

}
