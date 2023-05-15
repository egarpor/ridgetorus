
S <- function(s1, s2, rho) {
  rbind(c(s1^2, s1 * s2 * rho),
        c(s1 * s2 * rho, s2^2))
}
cmod <- function(x) (x + pi) %% (2 * pi) - pi
set.seed(1212)
x <- cmod(mvtnorm::rmvnorm(n = 200, sigma = S(1, 0.65, 0.9)))

pca <- princomp(x)
pc1 <- function(x) {
  m <- pca$loadings[2, 1] / pca$loadings[1, 1]
  cmod((pca$center[2] - m * pca$center[1]) + m * (x + pca$center[1]))
}

pdf("figures_wrapping/wrap1.pdf", width = 7, height = 7)
plot(x, xlim = c(-pi, pi), ylim = c(-pi, pi), axes = FALSE)
sdetorus::torusAxis()
y <- 1 * seq(-pi, pi, l = 1e5)
sdetorus::linesTorus(x = cmod(y), y = pc1(y), col = 2)
dev.off()

pdf("figures_wrapping/wrap2.pdf", width = 7, height = 7)
plot(x, xlim = c(-pi, pi), ylim = c(-pi, pi), axes = FALSE)
sdetorus::torusAxis()
y <- 2 * seq(-pi, pi, l = 1e5)
sdetorus::linesTorus(x = cmod(y), y = pc1(y), col = 2)
dev.off()

pdf("figures_wrapping/wrap4.pdf", width = 7, height = 7)
plot(x, xlim = c(-pi, pi), ylim = c(-pi, pi), axes = FALSE)
sdetorus::torusAxis()
y <- 4 * seq(-pi, pi, l = 1e5)
sdetorus::linesTorus(x = cmod(y), y = pc1(y), col = 2)
dev.off()

pdf("figures_wrapping/wrap100.pdf", width = 7, height = 7)
plot(x, xlim = c(-pi, pi), ylim = c(-pi, pi), axes = FALSE)
sdetorus::torusAxis()
y <- 100 * seq(-pi, pi, l = 1e5)
sdetorus::linesTorus(x = cmod(y), y = pc1(y), col = 2)
dev.off()

