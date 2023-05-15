library(ks)

# Gradient of a N(mu, Sigma) density (vectorized on x)
grad_norm <- function(x, mu, Sigma) {

  # Check dimensions
  x <- rbind(x)
  p <- length(mu)
  stopifnot(ncol(x) == p & nrow(Sigma) == p & ncol(Sigma) == p)

  # Gradient
  grad <- -mvtnorm::dmvnorm(x = x, mean = mu, sigma = Sigma) *
    t(t(x) - mu) %*% solve(Sigma)
  return(grad)

}

# Hessian of a N(mu, Sigma) density (vectorized on x)
Hess_norm <- function(x, mu, Sigma) {

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

  # As an array
  return(array(data = c(H), dim = c(p, p, nrow(x))))

}

# Projected gradient into the Hessian s-th eigenvector subspace
proj_grad_norm <- function(x, mu, Sigma, s = 2) {

  # Gradient
  grad <- grad_norm(x = x, mu = mu, Sigma = Sigma)

  # Hessian
  Hess <- Hess_norm(x = x, mu = mu, Sigma = Sigma)

  # Eigenvectors Hessian
  eig_Hess <- t(apply(Hess, 3, function(A) {
    eigen(x = A, symmetric = TRUE)$vectors[, s]
  }))

  # Projected gradient
  proj_grad <- t(sapply(1:nrow(eig_Hess), function(i) {
    tcrossprod(eig_Hess[i, ]) %*% grad[i, ]
  }))

  # As an array
  return(proj_grad)

}

# Compute vector fields
mu <- c(0, 0)
sigma_1 <- sqrt(1.25)
sigma_2 <- sqrt(2.5)
rho <- -0.5
Sigma <- matrix(c(sigma_1^2, sigma_1 * sigma_2 * rho,
                  sigma_1 * sigma_2 * rho, sigma_2^2), nrow = 2, ncol = 2)
x <- seq(-4, 4, l = 13)
xx <- as.matrix(expand.grid(x, x))
H <- Hess_norm(x = xx, Sigma = Sigma, mu = mu)
eig_val_1 <- t(apply(H, 3, function(A)
  eigen(x = A, symmetric = TRUE)$values[1]))
eig_val_2 <- t(apply(H, 3, function(A)
  eigen(x = A, symmetric = TRUE)$values[2]))
eig_vec_1 <- t(apply(H, 3, function(A)
  eigen(x = A, symmetric = TRUE)$vectors[, 1]))
eig_vec_2 <- t(apply(H, 3, function(A)
  eigen(x = A, symmetric = TRUE)$vectors[, 2]))
grad <- grad_norm(x = xx, mu = mu, Sigma = Sigma)
grad_proj <- proj_grad_norm(x = xx, mu = mu, Sigma = Sigma)

# Standardize directions for nicer plots
r <- 0.4
eig_vec_1 <- r * eig_vec_1
eig_vec_2 <- r * eig_vec_2
grad <- r * grad / sqrt(rowSums(grad^2))
grad_proj <- r * grad_proj / sqrt(rowSums(grad_proj^2))

pdf("figures_vectorfields/u1.pdf", width = 7, height = 7)
ks::plotmixt(mus = mu, Sigmas = Sigma, props = 1, display = "filled.contour2",
             gridsize = rep(251, 2), xlim = c(-4, 4), ylim = c(-4, 4),
             cont = seq(0, 95, by = 5), col.fun = viridis::viridis,
             binned = FALSE, xlab = "", ylab = "")
arrows(x0 = xx[, 1], y0 = xx[, 2],
       x1 = xx[, 1] + eig_vec_1[, 1], y1 = xx[, 2] + eig_vec_1[, 2],
       length = 0.1, angle = 20, col = 1, lwd = 3)
points(xx, pch = ifelse(eig_val_1 >= 0, "+", "-"),
       col = ifelse(eig_val_1 >= 0, 2, 4), cex = 2)
dev.off()

pdf("figures_vectorfields/u2.pdf", width = 7, height = 7)
ks::plotmixt(mus = mu, Sigmas = Sigma, props = 1, display = "filled.contour2",
             gridsize = rep(251, 2), xlim = c(-4, 4), ylim = c(-4, 4),
             cont = seq(0, 95, by = 5), col.fun = viridis::viridis,
             binned = FALSE, xlab = "", ylab = "")
arrows(x0 = xx[, 1], y0 = xx[, 2],
       x1 = xx[, 1] + eig_vec_2[, 1], y1 = xx[, 2] + eig_vec_2[, 2],
       length = 0.1, angle = 20, col = 1, lwd = 3)
points(xx, pch = ifelse(eig_val_2 >= 0, "+", "-"),
       col = ifelse(eig_val_2 >= 0, 2, 4), cex = 2)
dev.off()

pdf("figures_vectorfields/Df.pdf", width = 7, height = 7)
ks::plotmixt(mus = mu, Sigmas = Sigma, props = 1, display = "filled.contour2",
             gridsize = rep(251, 2), xlim = c(-4, 4), ylim = c(-4, 4),
             cont = seq(0, 95, by = 5), col.fun = viridis::viridis,
             binned = FALSE, xlab = "", ylab = "")
arrows(x0 = xx[, 1], y0 = xx[, 2],
       x1 = xx[, 1] + grad[, 1], y1 = xx[, 2] + grad[, 2],
       length = 0.1, angle = 20, col = 1, lwd = 3)
dev.off()

pdf("figures_vectorfields/Dfp.pdf", width = 7, height = 7)
ks::plotmixt(mus = mu, Sigmas = Sigma, props = 1, display = "filled.contour2",
             gridsize = rep(251, 2), xlim = c(-4, 4), ylim = c(-4, 4),
             cont = seq(0, 95, by = 5), col.fun = viridis::viridis,
             binned = FALSE, xlab = "", ylab = "")
arrows(x0 = xx[, 1], y0 = xx[, 2],
       x1 = xx[, 1] + grad_proj[, 1], y1 = xx[, 2] + grad_proj[, 2],
       length = 0.1, angle = 20, col = 1, lwd = 3)
dev.off()
