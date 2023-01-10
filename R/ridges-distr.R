
#' @noRd
#' @title Second eigenvalue of the Hessian
#'
#' @description Computation of the second eigenvalue of the Hessian for
#' bivariate von Mises and wrapped Cauchy densities.
#'
#' @param theta the points where to evaluate the density.
#' @inheritParams d_bwc
#' @inheritParams bvm
#' @inheritParams d_bwn
#' @param type either \code{"bvm"} for the bivariate von Mises sine or
#' \code{"bwc"} for the bivariate wrapped Cauchy
#' @return The value of the second eigenvalue at the point with the given
#' parameters.
#' @examples
#' xi1 <- 0.3
#' xi2 <- 0.5
#' rho <- 0.4
#' xi <- c(xi1, xi2, rho)
#' nth <- 50
#' th <- seq(-pi, pi, l = nth)
#' x <- as.matrix(expand.grid(th, th))
#' d <- ridgetorus:::H_eigenval_2(theta = x, xi = xi, type = "bwc")
#' filled.contour(th, th, matrix(d, nth, nth), col = viridisLite::viridis(20),
#'                levels = seq(min(d), max(d), l = 20))
H_eigenval_2 <- function(kmax = 1, theta, kappa, xi, mu, Sigma, type) {

  theta <- rbind(theta)

  if (type == "bvm") {

    k1 <- kappa[1]
    k2 <- kappa[2]
    lambda <- kappa[3]

    # Organize different terms
    sin1 <- sin(theta[, 1])
    cos1 <- cos(theta[, 1])
    sin2 <- sin(theta[, 2])
    cos2 <- cos(theta[, 2])
    sin12 <- sin1 * sin2

    # First derivatives
    C1 <- -k1 * sin1 + lambda * cos1 * sin2
    C2 <- -k2 * sin2 + lambda * sin1 * cos2
    C3 <- -k1 * cos1 - lambda * sin12
    C4 <- -k2 * cos2 - lambda * sin12
    C5 <- lambda * cos1 * cos2

    # Density value
    fvalue <- exp(k1 * cos1 + k2 * cos2 + lambda * sin12)

    # Second derivatives
    u <- (C3 + C1^2) * fvalue
    w <- (C4 + C2^2) * fvalue
    v <- (C5 + C1 * C2) * fvalue

    # Second eigenvalue
    rootd <- sqrt((w - u)^2 + 4 * v^2)
    return((u + w - rootd) / 2)

  } else if (type == "bwc") {

    xi1 <- xi[1]
    xi2 <- xi[2]
    rho <- xi[3]

    # Define constants
    rho2 <- rho^2
    xi12 <- xi1^2
    xi22 <- xi2^2
    c0 <- (1 + rho2) * (1 + xi12) * (1 + xi22) - 8 * abs(rho) * xi1 * xi2
    c1 <- 2 * (1 + rho2) * xi1 * (1 + xi22) - 4 * abs(rho) * (1 + xi12) * xi2
    c2 <- 2 * (1 + rho2) * (1 + xi12) * xi2 - 4 * abs(rho) * xi1 * (1 + xi22)
    c3 <- -4 * (1 + rho2) * xi1 * xi2 + 2 * abs(rho) * (1 + xi12) * (1 + xi22)
    c4 <- 2 * rho * (1 - xi12) * (1 - xi22)
    sin1 <- sin(-theta[, 1])
    cos1 <- cos(theta[, 1])
    sin2 <- sin(-theta[, 2])
    cos2 <- cos(theta[, 2])
    c1sin1 <- c1 * sin1
    c2sin2 <- c2 * sin2
    sin12 <- sin1 * sin2
    cos12 <- cos1 * cos2
    sin1cos2 <- sin1 * cos2
    sin2cos1 <- sin2 * cos1
    denominatorsq <- (c0 - c1 * cos1 - c2 * cos2 - c3 * cos12 - c4 * sin12)^2

    # Second derivatives
    u <- (c1sin1 + c3 * sin1cos2 - c4 * sin2cos1) *
      (2 * c1sin1 + 2 * c3 * sin1cos2 - 2 * c4 * sin2cos1) /
      (c0 - c1 * cos1 - c2 * cos2 - c3 * cos12 - c4 * sin12)^3 +
      (-c1 * cos1 - c3 * cos12 - c4 * sin12) / denominatorsq
    w <- (c2sin2 + c3 * sin2cos1 - c4 * sin1cos2) *
      (2 * c2sin2 + 2 * c3 * sin2cos1 - 2 * c4 * sin1cos2) /
      (c0 - c1 * cos1 - c2 * cos2 - c3 * cos12 - c4 * sin12)^3 +
      (-c2 * cos2 - c3 * cos12 - c4 * sin12) / denominatorsq
    v <- (c3 * sin12 + c4 * cos12) / denominatorsq +
      (c1sin1 + c3 * sin1cos2 - c4 * sin2cos1) *
      (2 * c2sin2 + 2 * c3 * sin2cos1 - 2 * c4 * sin1cos2) /
      (c0 - c1 * cos1 - c2 * cos2 - c3 * cos12 - c4 * sin12)^3

    # Second eigenvalue
    rootd <- sqrt((w - u)^2 + 4 * v^2)
    return((u + w - rootd) / 2)

  } else if (type == "bwn") {

    k <- seq.int(-kmax, kmax)
    K <- expand.grid(k, k)
    n <- nrow(theta)

    u <- numeric(length = n)
    v <- numeric(length = n)
    w <- numeric(length = n)
    Hess_norm <- function(x) {

      # Check dimensions
      x <- rbind(x)
      p <- length(mu)
      stopifnot(ncol(x) == p & nrow(Sigma) == p & ncol(Sigma) == p)

      # Hessian
      Sigma_inv <- solve(Sigma)
      H <- apply(x, 1, function(y) {
        mvtnorm::dmvnorm(x = y, mean = mu, sigma = Sigma) *
          (Sigma_inv %*% tcrossprod(y - mu) %*% Sigma_inv - Sigma_inv)
      })
      return(H)

    }

    for (i in seq_len(nrow(theta))) {

      th_expand <- 2 * pi * K + theta[i, ]
      Hessian <- Hess_norm(x = th_expand)
      u[i] <- sum(Hessian[1, ])
      v[i] <- sum(Hessian[2, ])
      w[i] <- sum(Hessian[4, ])

    }

    # Second eigenvalue
    rootd <- sqrt((w - u)^2 + 4 * v^2)
    return((u + w - rootd) / 2)

  }

}


#' @noRd
#' @title Ridge connected component that goes through the mode
#'
#' @description Computation of the connected component that goes through the
#' mode of the density ridge.
#'
#' @param solution the full set of points belonging to the ridge.
#' @param quadrant first or second quadrant.
#' @param max_dist maximum distance to be considered between points.
#' @return The points of the connected component of the ridge that goes through
#' the mode.
#' @examples
#' mu <- c(0, 0)
#' xi1 <- 0.3
#' xi2 <- 0.5
#' rho <- 0.4
#' nth <- 50
#' th <- seq(-pi, pi, l = nth)
#' x <- as.matrix(expand.grid(th, th))
#' d <- d_bwc(x = x, mu = mu, xi = c(xi1, xi2, rho))
#' filled.contour(th, th, matrix(d, nth, nth), col = viridisLite::viridis(20),
#'                levels = seq(0, max(d), l = 20))
filtering <- function(solution, quadrant, max_dist = 0.1) {

  # Take points in the given quadrant
  if (quadrant == 1) {

    solution <- solution[solution[, 1] >= 0, ]
    solution <- solution[solution[, 2] >= 0, ]

  }else if (quadrant == 2) {

    solution <- solution[solution[, 1] <= 0, ]
    solution <- solution[solution[, 2] >= 0, ]

  }

  theta1 <- solution[, 1]
  theta2 <- solution[, 2]

  # Initialize points in the connected component of the ridge
  theta1def <- c(0)
  theta2def <- c(0)
  lastpoint <- c(0, 0)

  # Avoid running out of points
  ntotal <- nrow(solution)
  j <- 1

  # While the closest point is near enough, keep adding points
  while (min(sqrt((theta1 - lastpoint[1])^2 + (theta2 - lastpoint[2])^2)) <
         max_dist && j < (ntotal - 1)) {

    # Calculate all distances
    distances <- sqrt((theta1 - lastpoint[1])^2 + (theta2 - lastpoint[2])^2)

    # Find the index of the point with the minimum distance
    index <- which.min(distances)

    # Find the closest point and update last_point
    closest_point <- solution[index, ]
    lastpoint <- closest_point

    # Update definitive list of points
    theta1def <- append(theta1def, closest_point[1])
    theta2def <- append(theta2def, closest_point[2])

    # Remove the point from the initial matrix
    theta1 <- theta1[-index]
    theta2 <- theta2[-index]
    solution <- solution[-index, ]
    j <- j + 1

  }

  # Calculate the remaining points by symmetry (excluding (0, 0))
  theta1def <- append(theta1def, -theta1def[-1])
  theta2def <- append(theta2def, -theta2def[-1])

  # Return a matrix with the connected component
  filtered_matrix <- matrix(nrow = length(theta1def), ncol = 2)
  filtered_matrix[, 1] <- theta1def
  filtered_matrix[, 2] <- theta2def
  return(filtered_matrix)

}


#' @title Connected component of the toroidal density ridge of a bivariate sine
#' von Mises, bivariate wrapped Cauchy, and bivariate wrapped normal
#'
#' @description Computation of the connected component of the density ridge of
#' in a given set of points or, if not specified, in a regular grid on
#' \eqn{[-\pi, \pi)}.
#'
#' @inheritParams bvm
#' @inheritParams bwc
#' @inheritParams bwn
#' @param eval_points evaluation points for the ridge.
#' @param subint_1 number of points for \eqn{\theta_1}.
#' @param subint_2 number of points for \eqn{\theta_2} at each \eqn{\theta_1}.
#' @return A matrix of size \code{c(subint_1, 2)} containing the points of the
#'  connected component of the ridge.
#' @references
#' Ozertem, U. and Erdogmus, D. (2011). Locally defined principal curves and
#' surfaces. \emph{Journal of Machine Learning Research}, 12(34):1249--1286.
#' \doi{10.6083/M4ZG6Q60}
#' @examples
#' \donttest{
#' # Bivariate von Mises
#' mu <- c(0, 0)
#' kappa <- c(0.3, 0.5, 0.4)
#' nth <- 100
#' th <- seq(-pi, pi, l = nth)
#' x <- as.matrix(expand.grid(th, th))
#' d <- d_bvm(x = x, mu = mu, kappa = kappa)
#' image(th, th, matrix(d, nth, nth), col = viridisLite::viridis(20))
#' ridge <- ridge_bvm(mu = mu, kappa = kappa, subint_1 = 5e2,
#'                    subint_2 = 5e2)
#' points(ridge)
#'
#' # Bivariate wrapped Cauchy
#' mu <- c(0, 0)
#' xi <- c(0.3, 0.6, 0.25)
#' nth <- 100
#' th <- seq(-pi, pi, l = nth)
#' x <- as.matrix(expand.grid(th, th))
#' d <- d_bwc(x = x, mu = mu, xi = xi)
#' image(th, th, matrix(d, nth, nth), col = viridisLite::viridis(20))
#' ridge <- ridge_bwc(mu = mu, xi = xi, subint_1 = 5e2, subint_2 = 5e2)
#' points(ridge)
#'
#' # Bivariate wrapped normal
#' mu <- c(0, 0)
#' Sigma <- matrix(c(10, 3, 3, 5), nrow = 2)
#' nth <- 100
#' th <- seq(-pi, pi, l = nth)
#' x <- as.matrix(expand.grid(th, th))
#' d <- d_bwn(x = x, mu = mu, Sigma = Sigma)
#' image(th, th, matrix(d, nth, nth), col = viridisLite::viridis(20))
#' ridge <- ridge_bwn(mu = mu, Sigma = Sigma, subint_1 = 5e2,
#'                    subint_2 = 5e2)
#' points(ridge)}
#' @export
#' @name ridge_distr


#' @rdname ridge_distr
#' @export
ridge_bvm <- function(mu, kappa, eval_points, subint_1, subint_2) {

  if (length(mu) != 2) {

    stop("mu must be a vector of length 2")

  }

  if (length(kappa) != 3) {

    stop("kappa must be a vector of length 3")

  }

  if (!all(kappa[1:2] >= 0)) {

    stop("Concentration parameters must be positive")

  }

  if (subint_1 <= 0 || subint_2 <= 0) {

    stop("Introduce a strictly positive number of subintervals")

  }

  # Search with the lowest kappa variable
  ord <- order(kappa[1:2])
  k1 <- kappa[c(ord == 1, FALSE)]
  k2 <- kappa[c(ord == 2, FALSE)]
  lambda <- kappa[3]
  mu1 <- mu[ord == 1]
  mu2 <- mu[ord == 2]

  # Limit cases with explicit solution
  if ((k1 / k2) == 1) {

    if (missing(eval_points)) {

      theta1 <- seq(-pi, pi, l = subint_1)

    } else {

      theta1 <- eval_points

    }

    sol_matrix <- matrix(nrow = length(theta1), ncol = 2)
    sol_matrix[, 1] <- sdetorus::toPiInt(theta1 + mu1)
    sol_matrix[, 2] <- sdetorus::toPiInt(sign(lambda) * theta1 + mu2)
    return(sol_matrix[, ord])

  } else {

    # Search through theta1 for the possible solution(s) of theta2
    if (missing(eval_points)) {

      theta1 <- seq(-pi, pi, l = subint_1)

    } else {

      theta1 <- eval_points

    }
    solution_theta1 <- c()
    solution_theta2 <- c()

    for (i in seq_len(subint_1)) {

      # Solve ridge implicit equation for candidate points
      sol <- rootSolve::uniroot.all(implicit_equation, theta1 = theta1[i],
                                    kappa = c(k1, k2),
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

          isridge[j] <- H_eigenval_2(theta = c(theta1[i], sol[j]),
                                     kappa = c(k1, k2, lambda),
                                     type = "bvm") <= 0

        }

        sol <- sol[isridge == 1]
        if (length(sol) > 0) {

          sol <- sol[abs(implicit_equation(theta2 = sol, theta1 = theta1[i],
                                           kappa = c(k1, k2),
                                           Lambda = matrix(c(0, lambda,
                                                             lambda, 0),
                                                           nrow = 2, ncol = 2),
                                           density = "bvm")) < 1e-4]

        }

        # Store the results
        theta1sol <- rep(theta1[i], length(sol))
        solution_theta1 <- append(solution_theta1, theta1sol)
        solution_theta2 <- append(solution_theta2, sol)

      }

    }

    # Matrix with results
    sol_matrix <- matrix(nrow = length(solution_theta1), ncol = 2)
    sol_matrix[, 1] <- solution_theta1
    sol_matrix[, 2] <- solution_theta2

    # Use filtering function to obtain the connected component that goes
    # through the mode
    quadrant <- ifelse(lambda < 0, 2, 1)
    sol_matrix_filter <- filtering(sol_matrix, quadrant = quadrant)
    sol_matrix_filter[, 1] <- sdetorus::toPiInt(sol_matrix_filter[, 1] + mu1)
    sol_matrix_filter[, 2] <- sdetorus::toPiInt(sol_matrix_filter[, 2] + mu2)

    # Final results
    return(sol_matrix_filter[, ord])

  }

}


#' @rdname ridge_distr
#' @export
ridge_bwc <- function(mu, xi, eval_points, subint_1, subint_2) {

  if (length(mu) != 2) {

    stop("mu must be a vector of length 2")

  }

  if (length(xi) != 3) {

    stop("xi must be a vector of length 3")

  }

  if (!all(xi[1:2] >= 0 & xi[1:2] < 1)) {

    stop("Concentration parameters must be inside [0, 1)")

  }

  if (abs(xi[3]) >= 1) {

    stop("Dependence parameter must be inside (-1, 1)")

  }

  if (subint_1 <= 0 || subint_2 <= 0) {

    stop("Introduce a strictly positive number of subintervals")

  }

  ord <- order(xi[1:2])
  xi1 <- xi[c(ord == 1, FALSE)]
  xi2 <- xi[c(ord == 2, FALSE)]
  rho <- xi[3]
  mu1 <- mu[ord == 1]
  mu2 <- mu[ord == 2]

  if (missing(eval_points)) {

    theta1 <- seq(-pi, pi, l = subint_1)

  } else {

    theta1 <- eval_points

  }

  # Limit cases with explicit solution
  if ((xi1 / xi2) == 1) {

    sol_matrix <- matrix(nrow = length(theta1), ncol = 2)
    sol_matrix[, 1] <- sdetorus::toPiInt(theta1 + mu1)
    sol_matrix[, 2] <- sdetorus::toPiInt(sign(rho) * theta1 + mu2)
    return(sol_matrix[, ord])

  } else {

    # Search through theta1 for the possible solution(s) of theta2
    solution_theta1 <- c()
    solution_theta2 <- c()

    for (i in seq_len(subint_1)) {

      # Solve ridge implicit equation for candidate points
      sol <- rootSolve::uniroot.all(implicit_equation, theta1 = theta1[i],
                                    xi = c(xi1, xi2, rho),
                                    interval = c(-pi, pi), n = subint_2,
                                    tol = 1e-8, density = "bwc")

      # Check for NAs (no solution at that theta1)
      sol <- na.omit(sol)
      lsol <- length(sol)

      if (lsol >= 1) {

        isridge <- numeric(length = lsol)
        for (j in seq_len(lsol)) {

          # Check the eigenvalue condition for each candidate
          isridge[j] <- H_eigenval_2(theta = c(theta1[i], sol[j]),
                                     xi = c(xi1, xi2, rho), type = "bwc") <= 0

        }

        sol <- sol[isridge == 1]
        if (length(sol) > 0) {

          sol <- sol[abs(implicit_equation(theta2 = sol, theta1 = theta1[i],
                                           xi = c(xi1, xi2, rho),
                                           density = "bwc")) < 1e-4]

        }

        # Store the results
        theta1sol <- rep(theta1[i], length(sol))
        solution_theta1 <- append(solution_theta1, theta1sol)
        solution_theta2 <- append(solution_theta2, sol)

      }

    }

    # Matrix with results
    sol_matrix <- matrix(nrow = length(solution_theta1), ncol = 2)
    sol_matrix[, 1] <- solution_theta1
    sol_matrix[, 2] <- solution_theta2

    # Use filtering function to obtain the connected component that goes
    # through the mode
    quadrant <- ifelse(rho < 0, 2, 1)
    sol_matrix_filtered <- filtering(sol_matrix, quadrant = quadrant)

    # sdetorus maps exact pi into -pi, which causes weird points
    sol_matrix_filtered[, 1] <- sdetorus::toInt(sol_matrix_filtered[, 1] +
                                                  mu1, a = -pi, b = pi + 1e-4)
    sol_matrix_filtered[, 2] <- sdetorus::toInt(sol_matrix_filtered[, 2] +
                                                  mu2, a = -pi, b = pi + 1e-4)

    # Final results
    return(sol_matrix_filtered[, ord])

  }

}


#' @rdname ridge_distr
#' @export
ridge_bwn <- function(mu, Sigma, kmax = 2, eval_points,
                      subint_1, subint_2) {

  if (!is.matrix(Sigma)) {

    stop("Sigma needs to be a 2x2 matrix")

  }

  if (!all(dim(Sigma) == 2)) {

    stop("Sigma needs to be a 2x2 matrix")

  }

  if (length(mu) != 2) {

    stop("mu must be a vector of length 2")

  }

  if (subint_1 <= 0 || subint_2 <= 0) {

    stop("Introduce a strictly positive number of subintervals")

  }

  if (missing(eval_points)) {

    theta1 <- seq(-pi, pi, l = subint_1)

  } else {

    theta1 <- eval_points

  }
  k <- seq.int(-kmax, kmax)
  K <- as.matrix(expand.grid(k, k))

  # Search through theta1 for the possible solution(s) of theta2
  solution_theta1 <- c()
  solution_theta2 <- c()

  for (i in seq_len(subint_1)) {

    # Solve ridge implicit equation for candidate points
    sol <- rootSolve::uniroot.all(implicit_equation, theta1 = theta1[i],
                                  Sigma = Sigma, mu = mu, k = K,
                                  interval = c(-pi, pi), n = subint_2,
                                  tol = 1e-8, density = "bwn")

    # Check for NAs (no solution at that theta1)
    sol <- na.omit(sol)
    lsol <- length(sol)

    if (lsol >= 1) {

      isridge <- numeric(length = lsol)
      for (j in seq_len(lsol)) {

        # Check the eigenvalue condition for each candidate
        isridge[j] <- H_eigenval_2(theta = c(theta1[i], sol[j]), kmax = kmax,
                                   Sigma = Sigma, mu = mu, type = "bwn") <= 0

      }

      sol <- sol[isridge == 1]
      if (length(sol) > 0) {

        sol <- sol[abs(implicit_equation(theta2 = sol, theta1 = theta1[i],
                                       mu = mu, Sigma = Sigma, k = K,
                                       density = "bwn")) < 1e-4]

      }

      # Store the results
      theta1sol <- rep(theta1[i], length(sol))
      solution_theta1 <- append(solution_theta1, theta1sol)
      solution_theta2 <- append(solution_theta2, sol)

    }

  }

  # Matrix with results
  sol_matrix <- matrix(nrow = length(solution_theta1), ncol = 2)
  sol_matrix[, 1] <- solution_theta1
  sol_matrix[, 2] <- solution_theta2

  # Use filtering function to obtain the connected component that goes
  # through the mode
  quadrant <- ifelse(Sigma[2, 1] < 0, 2, 1)
  sol_matrix_filtered <- filtering(solution = sol_matrix, quadrant = quadrant)
  sol_matrix_filtered[, 1] <- sdetorus::toPiInt(sol_matrix_filtered[, 1] +
                                                  mu[1])
  sol_matrix_filtered[, 2] <- sdetorus::toPiInt(sol_matrix_filtered[, 2] +
                                                  mu[2])

  # Final results
  return(sol_matrix_filtered)

}
