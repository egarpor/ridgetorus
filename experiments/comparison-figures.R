
# Required packages
library(ridgetorus)
library(latex2exp)
library(mixtools)

# Set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Function that plots the density ridge given the fit
show_ridge <- function(fit, col_data = 1, main = "", N = 5e2) {

  # Extract data
  x <- fit$data

  # Plot data
  pch <- 16
  plot(x[ord, ], xlim = c(-pi, pi), ylim = c(-pi, pi), axes = FALSE,
       col = col_data, xlab = expression(theta[1]), ylab = expression(theta[2]),
       pch = pch, cex.lab = 1.25)
  title(main = main)
  sdetorus::torusAxis(cex.axis = 1.25)

  # Plot ridge curve
  th_grid <- arclength_ridge_curve(mu = fit$mu_hat, coefs = fit$coefs_hat,
                                   ind_var = fit$ind_var, N = N)
  th_grid <- c(th_grid, th_grid[1])
  y <- ridge_curve(theta = th_grid, mu = fit$mu_hat, coefs = fit$coefs_hat,
                   ind_var = fit$ind_var)
  sdetorus::linesTorus(x = y[, 1], y = y[, 2], col = "black", lwd = 3)

}

## Example with concentrated wrapped normal

# Simulate data
set.seed(1)
n <- 250
var1 <- 0.2
var2 <- 0.8
rho <- 0.25
Sigma <- matrix(c(var1, sqrt(var1 * var2) * rho, sqrt(var1 * var2) * rho, var2),
                nrow = 2, ncol = 2)
data <- sdetorus::toPiInt(rmvnorm(n = n, mu = c(-pi, 0), sigma = Sigma))

# APCA
mu1 <- circular:::MeanCircularRad(data[, 1])
mu2 <- circular:::MeanCircularRad(data[, 2])
data_apca <- sdetorus::toPiInt(cbind(data[, 1] - mu1, data[, 2] - mu2))
X_pcs_wc <- prcomp(data_apca, scale. = FALSE, center = FALSE)
v <- X_pcs_wc$rotation[, 1]
slope <- v[2] / v[1]
intercept <- 0
pev <- cumsum(X_pcs_wc$sdev^2) / sum(X_pcs_wc$sdev^2)
pev
var_exp <- frechet_ss(x = X_pcs_wc$x)$var_exp
var_exp
# The PC1 is (x, y) = (mu1 + v[1] * t, mu2 + v[2] * t), i.e.,
# y = (mu2 - v2 / v1 * mu1) + (v2 / v1) * x
# However, in APCA, (mu1, mu2) are set to zero since the data is
# previously centered by the circular means

# TR-PCA
trpca <- ridge_pca(data_apca, type = "auto")
trpca$var_exp

# Order the scores and assign rainbow colors
scores_conc <- trpca$scores
ord <- order(scores_conc[, 1])
colors <- rainbow(n)

# Save ridge, PCA result and original data in one picture
pdf(file = "figures_comparison/bwn_1_fits.pdf")
show_ridge(trpca, col_data = colors)
abline(a = intercept, b = slope, col = "red", lwd = 3)
dev.off()

# Plot TR-PCA scores
pdf(file = "figures_comparison/bwn_1_scores_trpca.pdf")
plot(trpca$scores[ord, ], col = colors, pch = 16, xlim = c(-pi, pi),
     axes = FALSE, ylim = c(-pi, pi), xlab = TeX("$s_1$"), ylab = TeX("$s_2$"),
     cex.lab = 1.25)
sdetorus::torusAxis(cex.axis = 1.25)
abline(h = 0, lty = 3)
dev.off()

# Plot APCA scores
range(X_pcs_wc$x)
X_pcs_wc$x[, 1] <- -1 * X_pcs_wc$x[, 1]
pdf(file = "figures_comparison/bwn_1_scores_apca.pdf")
plot(X_pcs_wc$x[ord, ], col = colors, pch = 16, xlim = c(-pi, pi),
     ylim = c(-pi, pi), xlab = TeX("$s_1$"), ylab = TeX("$s_2$"), axes = FALSE,
     cex.lab = 1.25)
sdetorus::torusAxis(cex.axis = 1.25)
abline(h = 0, lty = 3)
dev.off()

## Example with a more spread wrapped normal

# Simulate data
set.seed(2)
n <- 250
var1 <- 3
var2 <- 1.5
rho <- 0.85
Sigma <- matrix(c(var1, sqrt(var1 * var2) * rho, sqrt(var1 * var2) * rho, var2),
                ncol = 2, nrow = 2)
data <- r_bwn(n = n, mu = c(pi / 2, 0), Sigma = Sigma)

# APCA
mu1 <- circular:::MeanCircularRad(data[, 1])
mu2 <- circular:::MeanCircularRad(data[, 2])
data_apca <- sdetorus::toPiInt(cbind(data[, 1] - mu1, data[, 2] - mu2))
X_pcs_wc <- prcomp(data_apca, scale. = FALSE, center = FALSE)
v <- X_pcs_wc$rotation[, 1]
slope <- v[2] / v[1]
intercept <- 0
pev <- cumsum(X_pcs_wc$sdev^2) / sum(X_pcs_wc$sdev^2)
pev
var_exp <- frechet_ss(x = X_pcs_wc$x)$var_exp
var_exp

# TR-PCA
trpca <- ridge_pca(data_apca, type = "auto")
trpca$var_exp

# Order the scores and assign rainbow colors
scores_conc <- trpca$scores
ord <- order(scores_conc[, 1])
colors <- rainbow(n)

# Save ridge, PCA result and original data in one picture
par(cex = 1.2)
pdf(file = "figures_comparison/bwn_2_fits.pdf")
show_ridge(trpca, col_data = colors)
abline(a = intercept, b = slope, col = "red", lwd = 3)
dev.off()

# Plot TR-PCA scores
pdf(file = "figures_comparison/bwn_2_scores_trpca.pdf")
plot(trpca$scores[ord, ], col = colors, pch = 16, xlim = c(-pi, pi),
     ylim = c(-pi, pi), xlab = TeX("$s_1$"), ylab = TeX("$s_2$"), axes = FALSE,
     cex.lab = 1.25)
sdetorus::torusAxis(cex.axis = 1.25)
abline(h = 0, lty = 3)
dev.off()

# Plot APCA scores
range(X_pcs_wc$x)
pdf(file = "figures_comparison/bwn_2_scores_apca.pdf")
plot(X_pcs_wc$x[ord, ], col = colors, pch = 16, xlim = c(-4.5, 4.5),
     ylim = c(-4.5, 4.5), xlab = TeX("$s_1$"), ylab = TeX("$s_2$"),
     axes = FALSE, cex.lab = 1.25)
sdetorus::torusAxis(cex.axis = 1.25)
abline(h = 0, lty = 3)
dev.off()

## Example with wrapped Cauchy

# Simulate data
set.seed(3)
n <- 250
data <- r_bwc(n = n, mu = c(1, 2), xi = c(0.5, 0.1, -0.75))

# APCA
mu1 <- circular:::MeanCircularRad(data[, 1])
mu2 <- circular:::MeanCircularRad(data[, 2])
data_apca <- sdetorus::toPiInt(cbind(data[, 1] - mu1,
                                     data[, 2] - mu2))
X_pcs_wc <- prcomp(data_apca, scale. = FALSE, center = FALSE)
v <- X_pcs_wc$rotation[, 1]
slope <- v[2] / v[1]
intercept <- 0
pev <- cumsum(X_pcs_wc$sdev^2) / sum(X_pcs_wc$sdev^2)
pev
var_exp <- frechet_ss(x = X_pcs_wc$x)$var_exp
var_exp

# TR-PCA
trpca <- ridge_pca(data_apca, type = "auto")
trpca$var_exp

# Order the scores and assign rainbow colors
scores_conc <- trpca$scores
ord <- order(scores_conc[, 1])
colors <- rainbow(n)

# Save ridge, PCA result and original data in one picture
par(cex = 1.2)
pdf(file = "figures_comparison/bwc_fits.pdf")
show_ridge(trpca, col_data = colors)
abline(a = intercept, b = slope, col = "red", lwd = 3)
dev.off()

# Plot TR-PCA scores
pdf(file = "figures_comparison/bwc_scores_trpca.pdf")
plot(trpca$scores[ord, ], col = colors, pch = 16, xlim = c(-pi, pi),
     ylim = c(-pi, pi), xlab = TeX("$s_1$"), ylab = TeX("$s_2$"), axes = FALSE,
     cex.lab = 1.25)
sdetorus::torusAxis(cex.axis = 1.25)
abline(h = 0, lty = 3)
dev.off()

# Plot APCA scores
range(X_pcs_wc$x)
X_pcs_wc$x[, 1] <- -1 * X_pcs_wc$x[, 1]
pdf(file = "figures_comparison/bwc_scores_apca.pdf")
plot(X_pcs_wc$x[ord, ], col = colors, pch = 16, xlim = c(-4.5, 4.5),
     ylim = c(-4.5, 4.5), xlab = TeX("$s_1$"), ylab = TeX("$s_2$"),
     cex.lab = 1.25, axes = FALSE)
sdetorus::torusAxis(cex.axis = 1.25)
abline(h = 0, lty = 3)
dev.off()

## Simpson's paradox example with two clusters

# Simulate data
set.seed(4)
n <- 250
var1 <- 0.4
var2 <- 0.16
rho <- 0.35
Sigma1 <- matrix(c(var1, sqrt(var1 * var2) * rho,
                   sqrt(var1 * var2) * rho, var2), nrow = 2, ncol = 2)
Sigma2 <- matrix(c(var2, sqrt(var1 * var2) * rho,
                   sqrt(var1 * var2) * rho, var1), nrow = 2, ncol = 2)
data <- rbind(r_bwn(n = n / 2, mu = c(pi / 2, -pi / 2), Sigma = Sigma1),
              r_bwn(n = n / 2, mu = c(-pi / 2, pi / 2), Sigma = Sigma2))

# APCA
mu1 <- circular:::MeanCircularRad(data[, 1])
mu2 <- circular:::MeanCircularRad(data[, 2])
data_apca <- sdetorus::toPiInt(cbind(data[, 1] - mu1, data[, 2] - mu2))
X_pcs_wc <- prcomp(data_apca, scale. = FALSE, center = FALSE)
v <- X_pcs_wc$rotation[, 1]
slope <- v[2] / v[1]
intercept <- 0
pev <- cumsum(X_pcs_wc$sdev^2) / sum(X_pcs_wc$sdev^2)
pev
var_exp <- frechet_ss(x = X_pcs_wc$x)$var_exp
var_exp

# TR-PCA
trpca <- ridge_pca(data_apca, type = "auto")
trpca$var_exp

# Order the scores and assign rainbow colors
scores_conc <- trpca$scores
ord <- order(scores_conc[, 1])
colors <- rep(3:4, each = n / 2)

# Save ridge, PCA result and original data in one picture
par(cex = 1.2)
pdf(file = "figures_comparison/simpson_fits.pdf")
show_ridge(trpca, col_data = colors)
abline(a = intercept, b = slope, col = "red", lwd = 3)
dev.off()

# Plot TR-PCA scores
pdf(file = "figures_comparison/simpson_scores_trpca.pdf")
plot(trpca$scores[ord, ], col = colors, pch = 16, xlim = c(-pi, pi),
     ylim = c(-pi, pi), xlab = TeX("$s_1$"), ylab = TeX("$s_2$"), axes = FALSE,
     cex.lab = 1.25)
sdetorus::torusAxis(cex.axis = 1.25)
abline(h = 0, lty = 3)
dev.off()

# Plot APCA scores
range(X_pcs_wc$x)
pdf(file = "figures_comparison/simpson_scores_apca.pdf")
plot(X_pcs_wc$x[ord, ], col = colors, pch = 16, xlim = c(-4.5, 4.5),
     ylim = c(-4.5, 4.5), xlab = TeX("$s_1$"), ylab = TeX("$s_2$"),
     axes = FALSE, cex.lab = 1.25)
sdetorus::torusAxis(cex.axis = 1.25)
abline(h = 0, lty = 3)
dev.off()
