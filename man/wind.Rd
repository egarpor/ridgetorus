% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{wind}
\alias{wind}
\title{Texas wind dataset}
\format{
A data frame with 30 rows and 2 variables:
\describe{
  \item{theta1}{Direction at 6:00 am.}
  \item{theta2}{Direction at 12:00 noon.}
}
}
\usage{
wind
}
\description{
Wind direction at 6:00 and 7:00 from June 1, 2003 to June 30,
2003, in radians, measured at a weather station in Texas coded as C28-1.
}
\details{
The direction is measured in radians in \eqn{[-\pi, \pi)} with
\eqn{-\pi}/\eqn{-\frac{-\pi}{2}}/\eqn{0}/\eqn{-\frac{\pi}{2}}/\eqn{\pi}
representing the East/South/West/North/East directions.
}
\examples{
# Load data
data("wind")
plot(wind, xlab = expression(theta[1]), ylab = expression(theta[2]),
     xlim = c(-pi, pi), ylim = c(-pi, pi), axes = FALSE)
sdetorus::torusAxis()

\donttest{
# Perform TR-PCA
fit <- ridge_pca(x = wind)
show_ridge_pca(fit)}
}
\references{
Johnson, R. A. and Wehrly, T. (1977). Measures and models for angular
correlation and angular-linear correlation. \emph{Journal of the Royal
Statistical Society. Series B (Methodological)}, 39(2):222--229.
\url{https://www.jstor.org/stable/2984799}
}
\keyword{datasets}
