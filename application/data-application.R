
# This script reproduces all the figures and analyses of the application of
# toroidal PCA via density ridges onto the Santa Barbara dataset.

# Required packages
library(ridgetorus)
library(ks)
library(latex2exp)
library(animation)
library(ggmap)
library(ggimage)
library(BBmisc)

# Set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load data
data("santabarbara")

# Custom filled.contour
my.filled.contour <- function(x = seq(0, 1, length.out = nrow(z)),
  y = seq(0, 1, length.out = ncol(z)), z, xlim = range(x, finite = TRUE),
  ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE),
  levels = pretty(zlim, nlevels), nlevels = 20,
  color.palette = function(n) hcl.colors(n, "YlOrRd", rev = TRUE),
  col = color.palette(length(levels) - 1), asp = NA, xaxs = "i", yaxs = "i",
  las = 1, axes = TRUE, frame.plot = axes, xlab = "x", ylab = "y", cexlab = 1,
  ...) {

  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      }
      else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    }
    else stop("no 'z' matrix specified")
  }
  else if (is.list(x)) {
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
  if (!missing(xlab) & !missing(ylab)){
    title(main = "", xlab = xlab, ylab = ylab, cex.lab = cexlab, ...)
  }

  .filled.contour(x, y, z, levels, col, ...)

}

# Function that given a 2D matrix with column names loc, performs the
# kernel density estimation, bvm and bwc fits, and ridge, as well as computes
# and plots the scores for the best fit
analysis <- function(x, loc, col = "red") {

  # Compute kde for a diagonal bandwidth matrix centering the data to mitigate
  # periodicity
  a_mean <- c(circular::mean.circular(x[, 1]), circular::mean.circular(x[, 2]))
  x_cent <- cbind(sdetorus::toPiInt(x[, 1] - a_mean[1]),
                  sdetorus::toPiInt(x[, 2] - a_mean[2]))
  Hns <- ks::Hns(x = x_cent)

  # Expand the grid to account for periodicity
  expanded_x <- 0
  for (i in c(-1, 0, 1)) {
    for (j in c(-1, 0, 1)) {

      expanded_x <- rbind(expanded_x, cbind(x[, 1] + i * 2 * pi ,
                                            x[, 2] + j * 2 * pi))

    }
  }
  expanded_x <- expanded_x[-1,]

  # Compute the kde
  kde <- ks::kde(x = expanded_x, H = Hns, gridsize = c(400, 400),
                 xmin = c(-3 * pi, -3 * pi), xmax = c(3 * pi, 3 * pi))

  ## BvM

  # Obtain the parameters of the von Mises distribution via ML
  fit <- ridge_pca(x = x, type = "bvm")
  param <- fit$fit_mle
  mu <- c(param$mu1, param$mu2)
  k <- c(param$kappa1, param$kappa2)
  lambda <- param$lambda

  # Print the parameters and the minus loglikelihood
  print(loc)
  print(fit$var_exp)
  print(mu)
  print(k)
  print(lambda)
  loglikelihood <- sum(log(d_bvm(x = x, mu = mu, kappa = c(k, lambda))))
  print(-loglikelihood)

  # Calculate density values and density ridge
  x1 <- seq( -pi, pi, l = 200)
  y1 <- seq( -pi, pi, l = 200)
  est <- apply(X = as.matrix(expand.grid(x1, y1)), 1, FUN = d_bvm,
               mu = mu, kappa = c(k, lambda))
  densityvalues_vm <- matrix(est, nrow = 200, ncol = 200, byrow = FALSE)
  ridge_vm <- ridge_bvm(mu = mu, kappa = c(k, lambda), subint_1 = 1000,
                        subint_2 = 1000)

  ## BWC

  # Obtain the parameters of the wrapped Cauchy distribution via ML
  fit <- ridge_pca(x = x, type = "bwc")
  param <- fit$fit_mle
  mu <- c(param$mu1, param$mu2)
  xi <- c(param$xi1, param$xi2)
  rho <- param$rho

  # Print the parameters and the minus log likelihood
  print(loc)
  print(fit$var_exp)
  print(mu)
  print(xi)
  print(rho)
  loglikelihood <- sum(log(d_bwc(x = x, mu = mu, xi = c(xi, rho))))
  print(-loglikelihood)

  # Calculate density values and density ridge
  est <- apply(X = as.matrix(expand.grid(x1, y1)), 1, FUN = d_bwc,
               mu = mu, xi = c(xi, rho))
  densityvalues_wc <- matrix(est, nrow = 200, ncol = 200, byrow = FALSE)
  ridge_wc <- ridge_bwc(mu = mu, xi = c(xi, rho), subint_1 = 1000,
                        subint_2 = 1000)

  ## Plotting

  # Choose the contour levels with the highest and lowest values of the density
  total_max <- max(c(max(densityvalues_wc), max(densityvalues_vm),
                     9 * max(kde$estimate)))
  total_min <- min(c(min(densityvalues_wc), min(densityvalues_vm),
                     9 * min(kde$estimate)))
  levels <- seq(total_min, total_max, l = 20)

  # Plot the results
  png(paste("figures_app/kde", loc[1], "_", loc[2], ".png", sep = ""),
      width = 7, height = 7, units = "in", res = 300, bg = "transparent")
  my.filled.contour(kde$eval.points[[1]], kde$eval.points[[2]],
                    9 * kde$estimate, col  =  viridis::viridis(23, alpha = 0.85),
                    axes = FALSE, xlab = TeX(paste("$\\theta_", loc[1],
                                                   sep = "")),
                    ylab = TeX(paste("$\\theta_", loc[2], sep = "")),
                    xlim = c(-pi, pi), ylim = c(-pi, pi), levels = levels,
                    cexlab = 1.25)
  at <- seq(-pi, pi, l = 5)
  labels <- c("W", "S", "E", "N", "W")
  axis(1, at = at, labels = labels, cex.axis = 1.25)
  axis(2, at = at, labels = labels, cex.axis = 1.25)
  points(x, cex = 0.5, pch = 16)
  dev.off()

  # Plot the results and save the image in a file
  png(paste("figures_app/vonMises", loc[1], "_", loc[2], ".png", sep = ""),
      width = 7, height = 7, units = "in", res = 300, bg = "transparent")
  my.filled.contour(x1, y1, densityvalues_vm, col = viridis::viridis(23, alpha = 0.85),
                    xlab = TeX(paste("$\\theta_", loc[1], sep = "")),
                    ylab = TeX(paste("$\\theta_", loc[2], sep = "")),
                    axes = FALSE, levels = levels, cexlab = 1.25)
  at <- seq(-pi, pi, l = 5)
  labels <- c("W", "S", "E", "N", "W")
  axis(1, at = at, labels = labels, cex.axis = 1.25)
  axis(2, at = at, labels = labels, cex.axis = 1.25)
  points(ridge_vm, col = col)
  dev.off()

  # Plot the results and save it into a file
  png(paste("figures_app/wrappedCauchy", loc[1], "_", loc[2], ".png", sep = ""),
      width = 7, height = 7, units = "in", res = 300, bg = "transparent")
  my.filled.contour(x1, y1, densityvalues_wc, col = viridis::viridis(23, alpha = 0.85),
                    xlab = TeX(paste("$\\theta_", loc[1], sep = "")),
                    ylab = TeX(paste("$\\theta_", loc[2], sep = "")),
                    axes = FALSE, levels = levels, cexlab = 1.25)
  at <- seq(-pi, pi, l = 5)
  labels <- c("W", "S", "E", "N", "W")
  axis(1, at = at, labels = labels, cex.axis = 1.25)
  axis(2, at = at, labels = labels, cex.axis = 1.25)
  points(ridge_wc, col = col)
  dev.off()

  # Plot final fit w/o ridge for the video
  fit <- ridge_pca(x = x, type = "auto")
  param <- fit$fit_mle

  if (fit$type == "bwc") {

    mu <- c(param$mu1, param$mu2)
    xi <- c(param$xi1, param$xi2)
    rho <- param$rho

    # Calculate density values and density ridge
    x1 <- seq(-pi, pi, l = 200)
    y1 <- seq(-pi, pi, l = 200)
    est <- apply(X = as.matrix(expand.grid(x1, y1)), 1, FUN = d_bwc,
                 mu = mu, xi = c(xi, rho))
    densityvalues_final <- matrix(est, nrow = 200, ncol = 200, byrow = FALSE)

  } else {

    mu <- c(param$mu1, param$mu2)
    k <- c(param$kappa1, param$kappa2)
    lambda <- param$lambda

    # Calculate density values and density ridge
    x1 <- seq( -pi, pi, l = 200)
    y1 <- seq( -pi, pi, l = 200)
    est <- apply(X = as.matrix(expand.grid(x1, y1)), 1, FUN = d_bvm,
                 mu = mu, kappa = c(k, lambda))
    densityvalues_final <- matrix(est, nrow = 200, ncol = 200, byrow = FALSE)

  }

  png(paste("figures_app/finalridge", loc[1], "_", loc[2], ".png", sep = ""),
      width = 7, height = 7, units = "in", res = 300, bg = "transparent")
  my.filled.contour(x1, y1, densityvalues_final, col = viridis::viridis(23),
                    axes = FALSE, levels = levels)
  dev.off()

  # Erase white spaces
  knitr::plot_crop(paste("figures_app/finalridge", loc[1],
                         "_", loc[2], ".png", sep = ""))

  ## Compute scores and their kde
  scores <- fit$scores

  # Compute kde for a diagonal bandwidth matrix
  Hns <- ks::Hns(x = scores)

  # Expand the grid to account for periodicity
  expanded_x <- 0
  for (i in c(-1, 0, 1)) {
    for (j in c(-1, 0, 1)) {

      expanded_x <- rbind(expanded_x, cbind(scores[, 1] + i * 2 * pi ,
                                            scores[, 2] + j * 2 * pi))

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

  png(paste("figures_app/scores", loc[1], "_", loc[2], ".png", sep = ""),
      width = 7, height = 7, units = "in", res = 300, bg = "transparent")
  my.filled.contour(kde$eval.points[[1]], kde$eval.points[[2]],
                    9 * kde$estimate, col  =  viridis::viridis(23, alpha = 0.85),
                    axes = FALSE, xlab = TeX("$\\s_1"),
                    ylab = TeX("$\\s_2"), levels = levels,
                    xlim = c(-pi, pi), ylim = c(-pi, pi), cexlab = 1.25)
  sdetorus::torusAxis(cex.axis = 1.25)
  points(scores, cex = 0.5, pch = 16)
  dev.off()

}

# Run the analysis
analysis(x = santabarbara[c("A", "B")], loc = c("A", "B"))
analysis(x = santabarbara[c("A", "C")], loc = c("A", "C"))
analysis(x = santabarbara[c("C", "D")], loc = c("C", "D"))

# Run the independence tests

# Implementation of phi^{(n)}(\lambda) based on Algorithm 1 in SM of
# "Nonparametric tests of independence for circular data based on
# trigonometric moments" (doi:10.5705/ss.202021.0416)
phi_lambda <- function(theta1, theta2, lambda = 0.1, B = 1e3, verbose = FALSE) {

  # Sample size
  n <- length(theta1)
  stopifnot(n == length(theta2))

  # As arrays for sphunif::Psi_mat
  dim(theta1) <- dim(theta2) <- c(n, 1, 1)

  # Prepare Jc1_mat and Jc1_mat computation
  ind_tri <- sphunif::upper_tri_ind(n = n)
  Jc <- function(theta) cos(lambda * sin(theta)) *
    exp(lambda * (cos(theta) - 1))
  # Jc <- function(theta, R = 1e3) {
  #   r <- 1:R
  #   2 * sum(cos(r * theta) / r^2)
  # }

  # Jc1_mat
  Jc1_mat <- matrix(0, nrow = n, ncol = n)
  upper_tri_R <- upper.tri(Jc1_mat)
  Jc1 <- Jc(sphunif::Psi_mat(data = theta1, ind_tri = ind_tri,
                             use_ind_tri = TRUE, angles_diff = TRUE))
  Jc1_mat[upper_tri_R] <- Jc1
  Jc1_mat <- Jc1_mat + t(Jc1_mat)
  diag(Jc1_mat) <- Jc(0)

  # Jc2_mat
  Jc2_mat <- matrix(0, nrow = n, ncol = n)
  Jc2 <- Jc(sphunif::Psi_mat(data = theta2, ind_tri = ind_tri,
                             use_ind_tri = TRUE, angles_diff = TRUE))
  Jc2_mat[upper_tri_R] <- Jc2
  Jc2_mat <- Jc2_mat + t(Jc2_mat)
  diag(Jc2_mat) <- Jc(0)

  # Original statistic. The diagonals of Jc1_mat and Jc2_mat are 1's
  const <- (1 / n^3) * (2 * sum(Jc1) + n) * (2 * sum(Jc2) + n)
  Tn <- (1 / n) * (2 * sum(Jc1 * Jc2) + n) + const -
    (2 / n^2) * sum(Jc1_mat %*% Jc2_mat)

  # Permutation loop
  if (verbose) pb <- txtProgressBar(style = 3)
  Tn_b <- numeric(B)
  for (b in 1:B) {

    # Randomly permute theta2
    theta2_b <- theta2[sample(n), , , drop = FALSE]

    # Jc2_mat_b
    Jc2_mat <- matrix(0, nrow = n, ncol = n)
    Jc2 <- Jc(sphunif::Psi_mat(data = theta2_b, ind_tri = ind_tri,
                               use_ind_tri = TRUE, angles_diff = TRUE))
    Jc2_mat[upper_tri_R] <- Jc2
    Jc2_mat <- Jc2_mat + t(Jc2_mat)
    diag(Jc2_mat) <- Jc(0)

    # Permuted statistic
    Tn_b[b] <- (1 / n) * (2 * sum(Jc1 * Jc2) + n) + const -
      (2 / n^2) * sum(Jc1_mat %*% Jc2_mat)

    # Progress
    if (verbose) setTxtProgressBar(pb, b / B)

  }

  # Test statistic
  Tn <- c("Tn" = Tn)

  # p-value
  p_value <- mean(Tn < Tn_b)

  # Test
  method <- "phi(lambda) test"
  test <- list(statistic = Tn, p.value = p_value,
               method = method, data.name = "(Theta_1, Theta_2)",
               boot_stats = Tn_b)
  class(test) <- "htest"
  return(test)

}

# Tests
phi_lambda(theta1 = santabarbara$A, theta2 = santabarbara$B, verbose = TRUE)
phi_lambda(theta1 = santabarbara$A, theta2 = santabarbara$C, verbose = TRUE)
phi_lambda(theta1 = santabarbara$C, theta2 = santabarbara$D, verbose = TRUE)

## Analysis for zone C-D

# Arc length of zone C-D

fit <- ridge_pca(x = santabarbara[c("C", "D")], type = "auto")
mu1 <- fit$mu_hat[1]
mu2 <- fit$mu_hat[2]

# Ridge is horizontal
order <- order(c(unname(unlist(fit$fit_mle[3])),
                 unname(unlist(fit$fit_mle[4]))))
coefs <- fit$coefs_hat

# Computes the arc length between t0 and t in intervals of length l
s_mod <- function(t, t0, l = 0.1) {

  # Avoid going backwards
  penalty <- 10

  dist <- dist_ridge_curve(alpha = c(t0, t), coefs = coefs)
  return(abs(dist - l) + penalty * (1 - sign(t - t0)))

}

# Solves s_mod to find equally spaced theta_1 in terms of arc length distance
arc_length <- function() {

  i <- 2
  t <- c(-pi)
  while (t[i - 1] < pi) {

    if (i > 3) {

      if (abs(t[i - 1] - t[i - 2]) < 0.01) {

        t <- append(t, optim(t[i - 1] + 0.015, s_mod, method = 'BFGS',
                             hessian = TRUE, t0 = t[i - 1])$par)

      } else {

        t <- append(t, optim(t[i - 1] + 0.005, s_mod, method = 'BFGS',
                             hessian = TRUE, t0 = t[i - 1])$par)

      }

    } else {

      t <- append(t, optim(t[i - 1] + 0.005, s_mod, method = 'BFGS',
                           hessian = TRUE, t0 = t[i - 1])$par)

    }

    i <- i + 1

  }

  return(t)

}

# Arcs
th1_arc <- arc_length()
th2_arc <- ridge_curve(theta = th1_arc, coefs = coefs)[, 2]

# Coordinates of the ridge in the miniplot. Image is 0.45 wide and 0.325 high

# Center of the square
mini_image_x <- -120.465
mini_image_y <- 34.55

# Size of the square
mini_image_width <- 0.45
mini_image_heigth <- 0.325

# Ridge image for the miniplot
d <- data.frame(x = mini_image_x, y = mini_image_y ,
                image = "figures_app/finalridgeC_D.png")

if (order[1] == 1) {

  theta_miniplot <- as.data.frame(cbind(sdetorus::toPiInt(th1_arc + mu1) *
                                          (mini_image_width + 0.009) / pi +
                                          mini_image_x - 0.011,
                                        sdetorus::toPiInt(th2_arc + mu2) *
                                          mini_image_heigth / pi +
                                          mini_image_y))

} else {

  theta_miniplot <- as.data.frame(cbind(sdetorus::toPiInt(th2_arc + mu1) *
                                          (mini_image_width + 0.009) / pi +
                                          mini_image_x - 0.011,
                                        sdetorus::toPiInt(th1_arc + mu2) *
                                          mini_image_heigth / pi +
                                          mini_image_y))

}

# Arcs
orderplot <- order(theta_miniplot[, 1])
th1_arc <- th1_arc[orderplot]
th2_arc <- th2_arc[orderplot]

# Arrow for the angle of the marine current at each zone A or zone B
arr_1 <- data.frame(u = cos(sdetorus::toPiInt(th1_arc + mu1)),
                    v = sin(sdetorus::toPiInt(th1_arc + mu1)))
arr_2 <- data.frame(u = cos(sdetorus::toPiInt(th2_arc + mu2)),
                    v = sin(sdetorus::toPiInt(th2_arc + mu2)))

# Load map of santabarbara and draw arrows and ridge on it
load("map_santa_barbara.RData")

# Coordinates of the two areas where the arrows should be drawn (C-D)
x_C <- -119.97
y_C <- 34.07
y_D <- 33.90
x_D <- -119.92
colors <- rainbow(n = length(th1_arc))

# Function to plot St. Barbara map, density ridge and current directions
ridge_index <- function(i, N_shadow = 10) {

  a <- ggmap(map.st_barbara) +
    theme(plot.margin = margin(0.1, 0.1, 0, 0, "cm"),
          panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent", colour = NA)) +
    geom_image(data = d, aes(x = x, y = y, image = image), size = .5) +
    geom_segment(aes(x = -120.904, y = 34.225, xend = -120.025,
                     yend = 34.225), col = "black") +
    geom_segment(aes(x = -120.025, y = 34.225, xend = -120.025,
                     yend = 34.87), col = "black") +
    geom_segment(aes(x = -120.904, y = 34.225, xend = -120.904,
                     yend = 34.87), col = "black") +
    geom_segment(aes(x = -120.904, y = 34.87, xend = -120.025,
                     yend = 34.87), col = "black") +
    xlab("Longitude") +
    ylab("Latitude") +
    # Draw ridge with rainbow
    geom_path(aes(x = V1, y = V2), data = theta_miniplot[orderplot, ],
              col = colors, linewidth = 1.15) +
    # Draw ball moving along the ridge
    geom_point(aes(x = V1, y = V2), data = theta_miniplot[orderplot, ][i, ],
               col = colors[i], size = 4)

  # Draw shadow
  if (N_shadow > 0) {

    # Transparency of the shadow arrows
    alpha <- seq(1, 0, l = N_shadow + 1)[-1] / 2

    # Draw arrows
    frame_seq <- seq( max(i - 1, 1), max(i - N_shadow, 1), by = -1)
    t <- 1
    for (j in frame_seq) {

      if (j != i) {

        a <- a + geom_segment(x = x_C, y = y_C, xend = x_C +
                                arr_1$u[j] / 8, yend = y_C + arr_1$v[j] / 8,
                              arrow = arrow(length = unit(0.45, "cm")),
                              col = colors[j], linewidth = 2, alpha = alpha[t])
        a <- a + geom_segment(x = x_D, y = y_D, xend = x_D +
                                arr_2$u[j] / 8, yend = y_D + arr_2$v[j] / 8,
                              arrow = arrow(length = unit(0.45, "cm")),
                              col = colors[j], linewidth = 2, alpha = alpha[t])
        t <- t + 1

      }

    }

    # Take into account periodicity for past arrows
    if (N_shadow > i) {

      for (j in seq(length(th1_arc), length(th1_arc) - N_shadow + i, by = -1)) {

        a <- a + geom_segment(x = x_C, y = y_C, xend = x_C +
                                arr_1$u[j] / 8, yend = y_C + arr_1$v[j] / 8,
                              arrow = arrow(length = unit(0.45, "cm")),
                              col = colors[j], linewidth = 2, alpha = alpha[t])
        a <- a + geom_segment(x = x_D, y = y_D, xend = x_D +
                                arr_2$u[j] / 8, yend = y_D + arr_2$v[j] / 8,
                              arrow = arrow(length = unit(0.45, "cm")),
                              col = colors[j], linewidth = 2, alpha = alpha[t])
        t <- t + 1

      }

    }

  }

  # Complete arrow on last layer on top
  a <- a +
    geom_segment(aes(x = x_C + 0.000001, y = y_C + 0.000001, xend = x_C +
                       arr_1$u[i] / 8, yend = y_C + arr_1$v[i] / 8),
                 arrow = arrow(length = unit(0.45, "cm")),
                 linewidth = 2.6, alpha = 1, linejoin = "round",
                 lineend = "round") +
    geom_segment(aes(x = x_C, y = y_C, xend = x_C +
                              arr_1$u[i] / 8, yend = y_C + arr_1$v[i] / 8),
                arrow = arrow(length = unit(0.45, "cm")), col = colors[i],
                linewidth = 2, alpha = 1, linejoin = "round",
                lineend = "round") +
    geom_segment(aes(x = x_D + 0.000001, y = y_D + 0.000001, xend = x_D +
                       arr_2$u[i] / 8, yend = y_D + arr_2$v[i] / 8),
                 arrow = arrow(length = unit(0.45, "cm")),
                 linewidth = 2.6, alpha = 1, linejoin = "round",
                 lineend = "round") +
    geom_segment(aes(x = x_D, y = y_D, xend = x_D +
                       arr_2$u[i] / 8, yend = y_D + arr_2$v[i] / 8),
                 arrow = arrow(length = unit(0.45, "cm")), col = colors[i],
                 linewidth = 2, alpha = 1, linejoin = "round",
                 lineend = "round")

  print(a)

}

# Set of 5 images
n <- length(th1_arc)
a <- 1
for (i in seq(5, n - 2, l = 5)) {

  print(a)
  png(paste("figures_app/indexing", a, ".png", sep = ""), width = 7, height = 7,
      units = "in", res = 150, bg = "transparent")
  ridge_index(i)
  dev.off()
  knitr::plot_crop(paste("figures_app/indexing", a, ".png", sep = ""))
  a <- a + 1

}

# Video
ani.options('interval' = 0.15, ani.width = 1400, ani.height = 1378,
            ani.res = 200, ani.dev = "png")
saveVideo(for (i in rep(seq(1, length(th1_arc), by = 1), 2)) ridge_index(i),
          video.name = "video_application.mp4")

