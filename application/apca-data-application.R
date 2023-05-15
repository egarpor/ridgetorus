
# Load libraries and data
library(ridgetorus)
library(dplyr)
library(latex2exp)
data("santabarbara")

## Plotting functions

my.filled.contour <- function(x = seq(0, 1, length.out = nrow(z)),
                              y = seq(0, 1, length.out = ncol(z)), z,
                              xlim = range(x, finite = TRUE),
                              ylim = range(y, finite = TRUE),
                              zlim = range(z, finite = TRUE),
                              levels = pretty(zlim, nlevels), nlevels = 20,
                              color.palette = function(n) hcl.colors(n, "YlOrRd",
                                                                     rev = TRUE),
                              col = color.palette(length(levels) - 1),
                              asp = NA, xaxs = "i", yaxs = "i", las = 1,
                              axes = TRUE, frame.plot = axes, xlab = "x",
                              ylab = "y", cexlab = 1,
                              ...) {

  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      } else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    } else {
      stop("no 'z' matrix specified")
    }
  } else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (any(diff(x) <= 0) || any(diff(y) <= 0))
    stop("increasing 'x' and 'y' values expected")
  mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
  on.exit(par(par.orig))
  w <- (3 + mar.orig[2L]) * par("csi") * 2.54

  mar <- mar.orig
  par(mar = mar)
  plot.new()
  plot.window(xlim, ylim, xaxs = xaxs, yaxs = yaxs, asp = asp
  )
  if (!missing(xlab) & !missing(ylab)) {
    title(main = "", xlab = xlab, ylab = ylab, cex.lab = cexlab, ...)
  }

  .filled.contour(x, y, z, levels, col, ...)

}

plot_scores <- function(scores, loc) {

  # Compute kde for a diagonal bandwidth matrix
  Hns <- ks::Hns(x = scores)

  # Compute the kde
  kde <- ks::kde(x = scores, H = Hns, gridsize = c(400, 400),
                 xmin = c(-3 * pi, -3 * pi), xmax = c(3 * pi, 3 * pi))

  # Solve ks bug that makes estimate 0
  minkde <- min(kde$estimate[kde$estimate > 1e-14])
  kde$estimate[kde$estimate < minkde] <- minkde
  levels <- seq(minkde, max(9 * kde$estimate), l = 20)

  png(paste("scoresapca", loc[1], "_", loc[2], ".png", sep = ""),
      width = 7, height = 7, units = "in", res = 300, bg = "transparent")
  my.filled.contour(kde$eval.points[[1]], kde$eval.points[[2]],
                    9 * kde$estimate,
                    col  =  viridis::viridis(23, alpha = 0.85),
                    axes = FALSE, xlab = TeX("$\\s_1$"),
                    ylab = TeX("$\\s_2$"), levels = levels,
                    xlim = c(-pi, pi) * 1.3, ylim = c(-pi, pi) * 1.3,
                    cexlab = 1.25)
  sdetorus::torusAxis(cex.axis = 1.25)
  points(scores, cex = 0.5, pch = 16)
  dev.off()

}

plot_apca <- function(data_apca, intercept, slope, loc) {

  # Compute kde for a diagonal bandwidth matrix
  Hns <- ks::Hns(x = data_apca)

  # Expand the grid to account for periodicity
  expanded_x <- 0
  for (i in c(-1, 0, 1)) {
    for (j in c(-1, 0, 1)) {

      expanded_x <- rbind(expanded_x, cbind(data_apca[, 1] + i * 2 * pi,
                                            data_apca[, 2] + j * 2 * pi))

    }
  }
  expanded_x <- expanded_x[-1, ]

  # Compute the kde
  kde <- ks::kde(x = expanded_x, H = Hns, gridsize = c(400, 400),
                 xmin = c(-3 * pi, -3 * pi), xmax = c(3 * pi, 3 * pi))

  # Solve ks bug that makes estimate 0
  minkde <- min(kde$estimate[kde$estimate > 1e-14])
  kde$estimate[kde$estimate < minkde] <- minkde
  levels <- seq(minkde, max(9 * kde$estimate), l = 20)

  # Compute TR-PCA results

  fit <- ridge_pca(x = data_apca, type = "auto")

  if (fit$type == "bwc") {

    param <- fit$fit_mle
    mu <- c(param$mu1, param$mu2)
    xi <- c(param$xi1, param$xi2)
    rho <- param$rho

    # Calculate density values and density ridge
    ridge <- ridge_bwc(mu = mu, xi = c(xi, rho), subint_1 = 1000,
                          subint_2 = 1000)

  } else {

    param <- fit$fit_mle
    mu <- c(param$mu1, param$mu2)
    k <- c(param$kappa1, param$kappa2)
    lambda <- param$lambda

    ridge <- ridge_bvm(mu = mu, kappa = c(k, lambda), subint_1 = 1000,
                          subint_2 = 1000)

  }

  png(paste("figures_aPCA/santabarbara_apca_", loc[1], "_", loc[2], ".png",
            sep = ""),
  width = 7, height = 7, units = "in", res = 300, bg = "transparent")
  my.filled.contour(kde$eval.points[[1]], kde$eval.points[[2]],
                    9 * kde$estimate,
                    col  =  viridis::viridis(23, alpha = 0.85),
                    axes = FALSE, xlab = TeX("$\\theta_1$"),
                    ylab = TeX("$\\theta_2$"), levels = levels,
                    xlim = c(-pi, pi), ylim = c(-pi, pi), cexlab = 1.25)
  sdetorus::torusAxis(cex.axis = 1.25)
  points(data_apca, cex = 0.5, pch = 16)
  abline(a = intercept, b = slope, col = "red", lwd = 3)
  points(ridge, col = "green", pch = 16, cex = 0.5)
  dev.off()

}

## A-B case

data <- cbind(santabarbara$A, santabarbara$B)
n <- length(data)
mu1 <- circular:::MeanCircularRad(data[, 1])
mu2 <- circular:::MeanCircularRad(data[, 2])
data_apca <- sdetorus::toPiInt(cbind(data[, 1] - mu1, data[, 2] - mu2))
X_pcs_wc <- prcomp(data_apca, scale. = FALSE, center = FALSE)
v <- X_pcs_wc$rotation[, 1]
slope <- v[2] / v[1]
intercept <- 0
var_exp <- frechet_ss(x = X_pcs_wc$x)$var_exp

# Plots both fits and scores aPCA
plot_apca(data_apca, intercept, slope, loc = c("A", "B"))
plot_scores(X_pcs_wc$x, loc = c("A", "B"))

# Order the scores and assign rainbow colors according to TR-PCA
trpca <- ridge_pca(data_apca, type = "auto")
scores_conc <- trpca$scores
ord <- order(scores_conc[, 1])
colors <- rainbow(dim(data_apca)[1])
png("figures_aPCA/rainbow_AB_TRPCA.png",
    width = 7, height = 7, units = "in", res = 300, bg = "transparent")
plot(data_apca[ord, ], cex = 0.5, pch = 16, col = colors, axes = FALSE,
     xlab = TeX("$\\theta_1$"), ylab = TeX("$\\theta_2$"))
sdetorus::torusAxis(cex.axis = 1.25)
dev.off()

# Order according to aPCA
ord <- order(X_pcs_wc$x[, 1])
colors <- rainbow(dim(data_apca)[1])
png("figures_aPCA/rainbow_AB_aPCA.png",
    width = 7, height = 7, units = "in", res = 300, bg = "transparent")
plot(data_apca[ord, ], cex = 0.5, pch = 16, col = colors, axes = FALSE,
     xlab = TeX("$\\theta_1$"), ylab = TeX("$\\theta_2$"))
sdetorus::torusAxis(cex.axis = 1.25)
dev.off()

## A-C case

data <- cbind(santabarbara$A, santabarbara$C)
mu1 <- circular:::MeanCircularRad(data[, 1])
mu2 <- circular:::MeanCircularRad(data[, 2])
data_apca <- sdetorus::toPiInt(cbind(data[, 1] - mu1, data[, 2] - mu2))
X_pcs_wc <- prcomp(data_apca, scale. = FALSE, center = FALSE)
v <- X_pcs_wc$rotation[, 1]
slope <- v[2] / v[1]
intercept <- 0
var_exp <- frechet_ss(x = X_pcs_wc$x)$var_exp

# Plots both fits and scores aPCA
plot_apca(data_apca, intercept, slope, loc = c("A", "C"))
plot_scores(X_pcs_wc$x, loc = c("A", "C"))

# Order the scores and assign rainbow colors according to TR-PCA
trpca <- ridge_pca(data_apca, type = "auto")
scores_conc <- trpca$scores
ord <- order(scores_conc[, 1])
colors <- rainbow(dim(data_apca)[1])
png("figures_aPCA/rainbow_AC_TRPCA.png",
    width = 7, height = 7, units = "in", res = 300, bg = "transparent")
plot(data_apca[ord, ], cex = 0.5, pch = 16, col = colors, axes = FALSE,
     xlab = TeX("$\\theta_1$"), ylab = TeX("$\\theta_2$"))
sdetorus::torusAxis(cex.axis = 1.25)
dev.off()

# Order according to aPCA
ord <- order(-X_pcs_wc$x[, 1])
colors <- rainbow(dim(data_apca)[1])
png("figures_aPCA/rainbow_AC_apca.png",
    width = 7, height = 7, units = "in", res = 300, bg = "transparent")
plot(data_apca[ord, ], cex = 0.5, pch = 16, col = colors, axes = FALSE,
     xlab = TeX("$\\theta_1$"), ylab = TeX("$\\theta_2$"))
sdetorus::torusAxis(cex.axis = 1.25)
dev.off()

## C-D case

data <- cbind(santabarbara$C, santabarbara$D)
mu1 <- circular:::MeanCircularRad(data[, 1])
mu2 <- circular:::MeanCircularRad(data[, 2])
data_apca <- sdetorus::toPiInt(cbind(data[, 1] - mu1, data[, 2] - mu2))
X_pcs_wc <- prcomp(data_apca, scale. = FALSE, center = FALSE)
v <- X_pcs_wc$rotation[, 1]
slope <- v[2] / v[1]
intercept <- 0
var_exp <- frechet_ss(x = X_pcs_wc$x)$var_exp

# Plots both fits and scores aPCA
plot_apca(data_apca, intercept, slope, loc = c("C", "D"))
plot_scores(X_pcs_wc$x, loc = c("C", "D"))

# Order the scores and assign rainbow colors according to TR-PCA
trpca <- ridge_pca(data_apca, type = "auto")
scores_conc <- trpca$scores
ord <- order(scores_conc[, 1])
colors <- rainbow(dim(data_apca)[1])

png("figures_aPCA/rainbow_CD_trpca.png",
    width = 7, height = 7, units = "in", res = 300, bg = "transparent")
plot(data_apca[ord, ], cex = 0.5, pch = 16, col = colors, axes = FALSE,
     xlab = TeX("$\\theta_1$"), ylab = TeX("$\\theta_2$"))
sdetorus::torusAxis(cex.axis = 1.25)
dev.off()

# Order according to aPCA
ord <- order(-X_pcs_wc$x[, 1])
colors <- rainbow(dim(data_apca)[1])
png("figures_aPCA/rainbow_CD_apca.png",
    width = 7, height = 7, units = "in", res = 300, bg = "transparent")
plot(data_apca[ord, ], cex = 0.5, pch = 16, col = colors, axes = FALSE,
     xlab = TeX("$\\theta_1$"), ylab = TeX("$\\theta_2$"))
sdetorus::torusAxis(cex.axis = 1.25)
dev.off()

## Performance comparison plots of TR-PCA vs aPCA

# Define area of interest
loc1 <- "A"; loc2 <- "B"
# loc1 <- "A"; loc2 <- "C"
# loc1 <- "C"; loc2 <- "D"

# Perform aPCA
data <- cbind(santabarbara[loc1], santabarbara[loc2])
n <- length(data)
mu1 <- circular:::MeanCircularRad(data[, 1])
mu2 <- circular:::MeanCircularRad(data[, 2])
data_apca <- sdetorus::toPiInt(cbind(data[, 1] - mu1, data[, 2] - mu2))
X_pcs_wc <- prcomp(data_apca, scale. = FALSE, center = FALSE)

# Define aPCA curve
v <- X_pcs_wc$rotation[, 1]
slope <- v[2] / v[1]
intercept <- 0

# TR-PCA
trpca <- ridge_pca(data_apca, type = "auto")

# Order according to the original distribution of points
ord <- order(data_apca[, 1])
colors <- rainbow(length(data_apca[, 1]))

# Plot the second scores vs toroidal distances between sample points and
# 1d-reconstructions of the sample points, rainbow-colored
# * For TR-PCA this is just scores_2 vs scores_2
# * For aPCA, this is the closest curve point to each sample point is needed

# Generate aPCA curve
apca_curve <- cbind(seq(-pi, pi, l = 1000), seq(-pi, pi, l = 1000) * slope)
apca_curve <- apca_curve[abs(apca_curve[, 2]) < pi, ]

# Find closest point in curve the each sample point via the minimum toroidal
# distance
ind_grid <- apply(data_apca, 1, function(th) {
  which.min(torus_dist(x = apca_curve, y = th, squared = TRUE))
})
theta_proj <- apca_curve[ind_grid, ]

# Toroidal distance to projection (unsigned)
pseudo_s2_apca_unsigned <- torus_dist(data_apca, theta_proj)

colors <- rainbow(nrow(X_pcs_wc$x))
ord <- order(data_apca[, 1])

png(paste("figures_aPCA/s2_vs_tdist_trpca_", loc1, "_", loc2, ".png", sep = ""),
    width = 7, height = 7, units = "in", res = 300, bg = "transparent")
plot(abs(trpca$scores[ord, 2]), abs(trpca$scores[ord, 2]), col = colors,
     axes = FALSE, xlab = TeX("$|\\s_2|$"),
     ylab = TeX("Toroidal distance to projection"),
     xlim = c(0, pi), ylim = c(0, pi))
sdetorus::torusAxis(cex.axis = 1.25)
dev.off()

png(paste("figures_aPCA/s2_vs_tdist_apca_", loc1, "_", loc2, ".png", sep = ""),
    width = 7, height = 7, units = "in", res = 300, bg = "transparent")
plot(abs(X_pcs_wc$x[ord, 2]), pseudo_s2_apca_unsigned[ord], col = colors,
     axes = FALSE, xlab = TeX("$|\\s_2|$"),
     ylab = TeX("Toroidal distance to projection"),
     xlim = c(0, pi), ylim = c(0, pi))
sdetorus::torusAxis(cex.axis = 1.25)
dev.off()
