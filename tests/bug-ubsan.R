library(ridgetorus)


# Bivariate wrapped normal
mu1 <- c(0, 0)
Sigma1 <- matrix(data = c(10, 2, 2, 5), nrow = 2)
debug(ridge_bwn)
ridge1 <- ridge_bwn(mu = mu1, Sigma = Sigma1, subint_1 = 5e2, subint_2 = 5e2)

ridgetorus:::implicit_equation(theta2 = numeric(0),
                               theta1 = -pi, mu = c(0, 0), Sigma = rbind(c(10, 2), c(2, 5)),
                               k = as.matrix(expand.grid(-2:2, -2:2)), density = "bwn")

