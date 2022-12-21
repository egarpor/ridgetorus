
# Required packages
library(ridgetorus)
library(latex2exp)

# Set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# BvM figures
bvm_plots <- function(mu = c(0, 0), kappa, subint_1, subint_2,
                      loc1 = TeX("$\\theta_1$"), loc2 = TeX("$\\theta_2$"),
                      col = "red", colors = viridis::viridis(22)) {

  ## Plot density on the background

  # x and y grid
  x1 <- seq(-pi, pi, l = 200)
  y1 <- seq(-pi, pi, l = 200)

  # Calculate density values
  est <- apply(X = as.matrix(expand.grid(x1, y1)), 1, FUN = d_bvm,
               kappa = kappa, mu = mu)
  densityvalues <- matrix(est, nrow = 200, ncol = 200, byrow = FALSE)

  # Plot density
  image(x1, y1, densityvalues, col = viridis::viridis(23), xlab = loc1,
        ylab = loc2, cex.lab = 1.25, axes = FALSE)
  sdetorus::torusAxis(cex.axis = 1.25)

  ## Solve in theta_2 for each theta_1 in seq(-pi, pi, l = 200)

  # Search with the lowest kappa variable
  ord <- order(kappa[1:2])
  k1 <- kappa[c(ord == 1, FALSE)]
  k2 <- kappa[c(ord == 2, FALSE)]
  lambda <- kappa[3]
  mu1 <- mu[ord == 1]
  mu2 <- mu[ord == 2]

  # Run solver in theta_2
  theta1 <- seq(-pi, pi, l = subint_1)
  solution_theta1 <- c()
  solution_theta2 <- c()
  for (i in seq_len(subint_1)) {

    # Solve ridge implicit equation for candidate points
    sol <- rootSolve::uniroot.all(ridgetorus:::implicit_equation,
                                  theta1 = theta1[i], kappa = c(k1, k2),
                                  Lambda = matrix(c(0, lambda, lambda, 0),
                                                  nrow = 2, ncol = 2),
                                  interval = c(-pi, pi), n = subint_2,
                                  tol = 1e-8, density = "bvm")

    # Check for NAs (no solution at that theta1)
    sol <- na.omit(sol)
    lsol <- length(sol)

    # Check the eigenvalue condition for each candidate
    if (lsol >= 1) {

      isridge <- numeric(length = lsol)

      for (j in seq_len(lsol)) {

        isridge[j] <- ridgetorus:::H_eigenval_2(theta = c(theta1[i], sol[j]),
                                                kappa = c(k1, k2, lambda),
                                                type = "bvm") <= 0

      }

      sol <- sol[isridge == 1]
      sol <- sol[abs(ridgetorus:::implicit_equation(
        theta2 = sol, theta1 = theta1[i], kappa = c(k1, k2),
        Lambda = matrix(c(0, lambda, lambda, 0), nrow = 2, ncol = 2),
        density = "bvm")) < 1e-4]

      # Store the results
      theta1sol <- rep(theta1[i], length(sol))
      solution_theta1 <- append(solution_theta1, theta1sol)
      solution_theta2 <- append(solution_theta2, sol)

    }
  }

  ## Solve in a perpendicular direction for higher resolution

  for (i in seq_len(subint_1)) {

    # Solve ridge implicit equation for candidate points
    sol <- rootSolve::uniroot.all(ridgetorus:::implicit_equation,
                                  theta1 = theta1[i], kappa = c(k2, k1),
                                  Lambda = matrix(c(0, lambda, lambda, 0),
                                                  nrow = 2, ncol = 2),
                                  interval = c(-pi, pi), n = subint_2,
                                  tol = 1e-8, density = "bvm")

    # Check for NAs (no solution at that theta1)
    sol <- na.omit(sol)
    lsol <- length(sol)

    # Check the eigenvalue condition for each candidate
    if (lsol >= 1) {

      isridge <- numeric(length = lsol)

      for (j in seq_len(lsol)) {

        isridge[j] <- ridgetorus:::H_eigenval_2(theta = c(theta1[i], sol[j]),
                                                kappa = c(k2, k1, lambda),
                                                type = "bvm") <= 0

      }

      sol <- sol[isridge == 1]
      sol <- sol[abs(ridgetorus:::implicit_equation(
        theta2 = sol, theta1 = theta1[i], kappa = c(k2, k1),
        Lambda = matrix(c(0, lambda, lambda, 0), nrow = 2, ncol = 2),
        density = "bvm")) < 1e-4]

      # Store the results
      theta1sol <- rep(theta1[i], length(sol))
      solution_theta1 <- append(solution_theta1, sol)
      solution_theta2 <- append(solution_theta2, theta1sol)

    }
  }

  ## Plot ridges

  # Plot ridge points
  sol_matrix <- matrix(nrow = length(solution_theta1), ncol = 2)
  sol_matrix[, 1] <- solution_theta1
  sol_matrix[, 2] <- solution_theta2
  points(sol_matrix[, ord], pch = 16)

  # Plot connected component
  if (k1 == k2) {

    sol_matrix_filter <- matrix(nrow = length(theta1), ncol = 2)
    sol_matrix_filter[, 1] <- theta1
    sol_matrix_filter[, 2] <- sign(lambda) * theta1

  } else {

    # Use filtering function to obtain the connected component that goes
    # through the mode
    quadrant <- ifelse(lambda < 0, 2, 1)
    sol_matrix_filter <- ridgetorus:::filtering(sol_matrix, quadrant = quadrant)
    sol_matrix_filter[, 1] <- sdetorus::toPiInt(sol_matrix_filter[, 1] + mu1)
    sol_matrix_filter[, 2] <- sdetorus::toPiInt(sol_matrix_filter[, 2] + mu2)

  }
  points(sol_matrix_filter[, ord], pch = 16, col = "red")

}

# BWC figures
bwc_plots <- function(mu = c(0, 0), xi, subint_1, subint_2,
                      loc1 = TeX("$\\theta_1$"), loc2 = TeX("$\\theta_2$"),
                      col = "red", colors = viridis::viridis(22)) {

  ## Plot density on the background

  # x and y grid
  x1 <- seq(-pi, pi, l = 200)
  y1 <- seq(-pi, pi, l = 200)

  # Calculate density values
  est <- apply(X = as.matrix(expand.grid(x1, y1)), 1, FUN = d_bwc, mu = mu,
               xi = xi)
  densityvalues <- matrix(est, nrow = 200, ncol = 200, byrow = FALSE)

  # Plot density
  image(x1, y1, densityvalues, col = viridis::viridis(23), xlab = loc1,
        ylab = loc2, cex.lab = 1.25, axes = FALSE)
  sdetorus::torusAxis(cex.axis = 1.25)

  ## Solve in theta_2 for each theta_1 in seq(-pi, pi, l = 200)

  # Search with the lowest kappa variable
  ord <- order(xi[1:2])
  xi1 <- xi[c(ord == 1, FALSE)]
  xi2 <- xi[c(ord == 2, FALSE)]
  rho <- xi[3]
  mu1 <- mu[ord == 1]
  mu2 <- mu[ord == 2]

  # Run solver in theta_2
  theta1 <- seq(-pi, pi, l = subint_1)
  solution_theta1 <- c()
  solution_theta2 <- c()
  for (i in seq_len(subint_1)) {

    # Solve ridge implicit equation for candidate points
    sol <- rootSolve::uniroot.all(ridgetorus:::implicit_equation,
                                  theta1 = theta1[i], xi = c(xi1, xi2, rho),
                                  interval = c(-pi, pi), n = subint_2,
                                  tol = 1e-8, density = "bwc")

    # Check for NAs (no solution at that theta1)
    sol <- na.omit(sol)
    lsol <- length(sol)

    if (lsol >= 1) {

      isridge <- numeric(length = lsol)

      for (j in seq_len(lsol)) {

        # Check the eigenvalue condition for each candidate
        isridge[j] <- ridgetorus:::H_eigenval_2(theta = c(theta1[i], sol[j]),
                                                xi = c(xi1, xi2, rho),
                                                type = "bwc") <= 0

      }

      sol <- sol[isridge == 1]
      sol <- sol[abs(ridgetorus:::implicit_equation(
        theta2 = sol, theta1 = theta1[i], xi = c(xi1, xi2, rho),
        density = "bwc")) < 1e-6]

      # Store the results
      theta1sol <- rep(theta1[i], length(sol))
      solution_theta1 <- append(solution_theta1, theta1sol)
      solution_theta2 <- append(solution_theta2, sol)

    }
  }

  ## Solve in a perpendicular direction for higher resolution

  for (i in seq_len(subint_1)) {

    # Solve ridge implicit equation for candidate points
    sol <- rootSolve::uniroot.all(ridgetorus:::implicit_equation,
                                  theta1 = theta1[i], xi = c(xi2, xi1, rho),
                                  interval = c(-pi, pi), n = subint_2,
                                  tol = 1e-8, density = "bwc")

    # Check for NAs (no solution at that theta1)
    sol <- na.omit(sol)
    lsol <- length(sol)

    if (lsol >= 1) {

      isridge <- numeric(length = lsol)

      for (j in seq_len(lsol)) {

        # Check the eigenvalue condition for each candidate
        isridge[j] <- ridgetorus:::H_eigenval_2(
          theta = c(theta1[i], sol[j]), xi = c(xi2, xi1, rho),
          type = "bwc") <= 0

      }

      sol <- sol[isridge == 1]
      sol <- sol[abs(ridgetorus:::implicit_equation(
        theta2 = sol, theta1 = theta1[i], xi = c(xi2, xi1, rho),
        density = "bwc")) < 1e-6]

      # Store the results
      theta1sol <- rep(theta1[i], length(sol))
      solution_theta1 <- append(solution_theta1, sol)
      solution_theta2 <- append(solution_theta2, theta1sol)

    }
  }

  ## Plot ridges

  # Plot ridge points
  solution_matrix <- matrix(nrow = length(solution_theta1), ncol = 2)
  solution_matrix[, 1] <- solution_theta1
  solution_matrix[, 2] <- solution_theta2
  points(solution_matrix[, ord], pch = 16)

  # Plot connected component
  if (xi1 == xi2) {

    sol_matrix_filter <- matrix(nrow = length(theta1), ncol = 2)
    sol_matrix_filter[, 1] <- theta1
    sol_matrix_filter[, 2] <- sign(rho) * theta1

  } else {

    # Use filtering function to obtain the connected component that goes
    # through the mode
    quadrant <- ifelse(sign(rho) == -1, 2, 1)
    sol_matrix_filter <- ridgetorus:::filtering(solution_matrix,
                                                quadrant = quadrant)

  }
  sol_matrix_filter[, 1] <- sdetorus::toPiInt(sol_matrix_filter[, 1] + mu1)
  sol_matrix_filter[, 2] <- sdetorus::toPiInt(sol_matrix_filter[, 2] + mu2)
  points(sol_matrix_filter[, ord], col = col, pch = 16)

}

## BvM plots

par(cex = 1.25)
png("figures_ridge/vonmises_03_015_025.png", width = 7, height = 7,
    units = "in", res = 300, bg = "transparent")
bvm_plots(kappa = c(0.3, 0.15, 0.25), subint_1 = 1000, subint_2 = 1000,
          colors = viridis::viridis(27))
dev.off()

png("figures_ridge/vonmises_03_06_05.png", width = 7, height = 7,
    units = "in", res = 300, bg = "transparent")
bvm_plots(kappa = c(0.3, 0.6, 0.5), subint_1 = 1000, subint_2 = 1000,
          colors = viridis::viridis(27))
dev.off()

png("figures_ridge/vonmises_03_03_1.png", width = 7, height = 7,
    units = "in", res = 300, bg = "transparent")
bvm_plots(kappa = c(0.3, 0.3, 1), subint_1 = 1000, subint_2 = 1000,
          colors = viridis::viridis(27))
dev.off()

png("figures_ridge/vonmises_1_05_15.png", width = 7, height = 7,
    units = "in", res = 300, bg = "transparent")
bvm_plots(kappa = c(1, 0.5, 1.5), subint_1 = 1000, subint_2 = 1000,
          colors = viridis::viridis(27))
dev.off()

## BWC plots

png("figures_ridge/bwc_horizontal.png", width = 7, height = 7,
    units = "in", res = 300, bg = "transparent")
bwc_plots(xi = c(0.2, 0.7, 0.2), subint_1 = 1000, subint_2 = 1000)
dev.off()

png("figures_ridge/bwc_Sdiagonal.png", width = 7, height = 7,
    units = "in", res = 300, bg = "transparent")
bwc_plots(xi = c(0.15, 0.075, 0.25), subint_1 = 1000, subint_2 = 1000)
dev.off()

png("figures_ridge/bwc_smallxi.png", width = 7, height = 7,
    units = "in", res = 300, bg = "transparent")
bwc_plots(xi = c(0.025, 0.6, 0.7), subint_1 = 1000, subint_2 = 1000)
dev.off()

png("figures_ridge/bwc_diagonal.png", width = 7, height = 7,
    units = "in", res = 300, bg = "transparent")
bwc_plots(xi = c(0.3, 0.3, 0.6), subint_1 = 1000, subint_2 = 1000)
dev.off()
