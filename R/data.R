

#' @title Santa Barbara currents
#'
#' @description The Santa Barbara Channel is a coastal area in California. This
#' dataset contains the sea currents in the four areas present in the
#' data application in García-Portugués and Prieto-Tirado (2022). Precisely,
#' it contains the 24-hour speed-weighted mean of the currents' direction in
#' each of the four areas downloaded from the
#' \href{https://hfradar.ndbc.noaa.gov/}{
#' NOAA High Frequency Radar National Server}.
#'
#' @docType data
#' @format A data frame with 1092 rows and 4 variables:
#' \describe{
#'   \item{A}{Sea current direction at zone A.}
#'   \item{B}{Sea current direction at zone B.}
#'   \item{C}{Sea current direction at zone C.}
#'   \item{D}{Sea current direction at zone D.}
#' }
#' @details
#' The selection of these four areas is motivated by previous studies on the
#' Santa Barbara currents, like Auad et al. (1998). The direction is measured
#' in radians in \eqn{[-\pi, \pi)} with \eqn{-\pi} / \eqn{-\frac{\pi}{2}} /
#' \eqn{0} / \eqn{\frac{\pi}{2}} / \eqn{\pi} representing the
#' East / South / West / North / East directions. The script performing the data
#' preprocessing is available at
#' \href{https://github.com/egarpor/ridgetorus/blob/master/data-raw/data-acquisition.R}{
#' \code{data-acquisition.R}}. The data was retrieved on 2022-10-21.
#' @references
#' Auad, G., Hendershott, M. C., and Winant, C. D. (1998). Wind-induced currents
#' and bottom-trapped waves in the Santa Barbara Channel. \emph{Journal of
#' Physical Oceanography}, 28(1):85--102.
#' \doi{10.1175/1520-0485(1998)028<0085:WICABT>2.0.CO;2}
#'
#' García-Portugués, E. and Prieto-Tirado, A. (2022). Toroidal PCA via density
#' ridges. \emph{arXiv:2212.10856}. \url{https://arxiv.org/abs/2212.10856}
#' @examples
#' # Load data
#' data("santabarbara")
#' AB_zone <- santabarbara[c("A","B")]
#'
#' \donttest{
#' # Perform TR-PCA
#' fit <- ridge_pca(x = AB_zone)
#' show_ridge_pca(fit)}
"santabarbara"


#' @title Japanese earthquakes dataset
#'
#' @description Pre-earthquake direction of steepest descent and the direction
#' of lateral ground movement before and after, respectively, an earthquake in
#' Noshiro (Japan) in 1983.
#'
#' @docType data
#' @format A data frame with 678 rows and 2 variables:
#' \describe{
#'   \item{theta1}{Direction of steepest descent.}
#'   \item{theta2}{Direction of lateral ground movement.}
#' }
#' @details
#' The direction is measured in radians in \eqn{[0, 2\pi)} with \eqn{0} /
#' \eqn{\frac{\pi}{2}} / \eqn{\pi} / \eqn{\frac{3\pi}{2}} / \eqn{2\pi}
#' representing the East / North / West / South / East directions.
#' @references
#' Hamada, M. and O'Rourke, T. (1992). Case Studies of Liquefaction & Lifeline
#' Performance During Past Earthquake. Volume 1: Japanese Case Studies.
#' Technical Report NCEER-92-0001. National Center for Earthquake Engineering
#' Research, University at Buffalo. \doi{10.1016/0886-7798(93)90146-m}
#'
#' Jones, M. C., Pewsey, A., and Kato, S. (2015). On a class of circulas:
#' copulas for circular distributions. \emph{Annals of the Institute of
#' Statistical Mathematics}, 67(5):843--862. \doi{10.1007/s10463-014-0493-6}
#'
#' Rivest, L.-P. (1997). A decentred predictor for circular-circular
#' regression. \emph{Biometrika}, 84(3):717--726. \doi{10.1093/biomet/84.3.717}
#' @examples
#' # Load data
#' data("earthquakes")
#'
#' # Transform the data into [-pi, pi)
#' earthquakes <- sdetorus::toPiInt(earthquakes)
#' plot(earthquakes, xlab = expression(theta[1]), ylab = expression(theta[2]),
#'      xlim = c(-pi, pi), ylim = c(-pi, pi), axes = FALSE)
#' sdetorus::torusAxis()
#'
#' \donttest{
#' # Perform TR-PCA
#' fit <- ridge_pca(x = earthquakes)
#' show_ridge_pca(fit)}
"earthquakes"


#' @title Texas wind dataset
#'
#' @description Wind direction at 6:00 and 7:00 from June 1, 2003 to June 30,
#' 2003, in radians, measured at a weather station in Texas coded as C28-1.
#'
#' @docType data
#' @format A data frame with 30 rows and 2 variables:
#' \describe{
#'   \item{theta1}{Direction at 6:00 am.}
#'   \item{theta2}{Direction at 12:00 noon.}
#' }
#' @details
#' The direction is measured in radians in \eqn{[-\pi, \pi)} with
#' \eqn{-\pi}/\eqn{-\frac{-\pi}{2}}/\eqn{0}/\eqn{-\frac{\pi}{2}}/\eqn{\pi}
#' representing the East/South/West/North/East directions.
#' @references
#' Johnson, R. A. and Wehrly, T. (1977). Measures and models for angular
#' correlation and angular-linear correlation. \emph{Journal of the Royal
#' Statistical Society. Series B (Methodological)}, 39(2):222--229.
#' \url{https://www.jstor.org/stable/2984799}
#' @examples
#' # Load data
#' data("wind")
#' plot(wind, xlab = expression(theta[1]), ylab = expression(theta[2]),
#'      xlim = c(-pi, pi), ylim = c(-pi, pi), axes = FALSE)
#' sdetorus::torusAxis()
#'
#' \donttest{
#' # Perform TR-PCA
#' fit <- ridge_pca(x = wind)
#' show_ridge_pca(fit)}
"wind"
