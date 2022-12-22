

#' @title \code{ridgetorus}: PCA on the Torus via Density Ridges
#'
#' @description Implementation of a Principal Component Analysis (PCA) in the
#' torus via density ridge estimation. The main function, ridge_pca(), obtains
#' the relevant density ridge for bivariate sine von Mises and bivariate wrapped
#' Cauchy distribution models and provides the associated scores and variance
#' decomposition. Auxiliary functions for evaluating, fitting, and sampling
#' these models are also provided. The package provides replicability to
#' García-Portugués and Prieto-Tirado (2022) <arXiv:2212.10856>.
#'
#' @author Eduardo García-Portugués and Arturo Prieto-Tirado.
#' @references
#' García-Portugués, E. and Prieto-Tirado, A. (2022). Toroidal PCA via density
#' ridges. \emph{arXiv:2212.10856}. \url{https://arxiv.org/abs/2212.10856}
#' @docType package
#' @name ridgetorus-package
#' @import graphics grDevices stats Rcpp
#' @useDynLib ridgetorus
#' @aliases ridgetorus ridgetorus-package
NULL
