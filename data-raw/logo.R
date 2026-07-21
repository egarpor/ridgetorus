
#
# Hex sticker logo for the ridgetorus package.
#
# The interior is the package's signature object: a bivariate wrapped Cauchy
# (BWC) density on the torus drawn as a viridis heatmap, with its density ridge
# overlaid in red. Mild and *unequal* marginal concentrations make the ridge a
# curve (rather than a straight diagonal); mild positive correlation tilts it.
#
# Re-run from the package root with:  Rscript data-raw/logo.R
# Output:                             man/figures/logo.png
#

library(ridgetorus)
library(ggplot2)
library(hexSticker)

## BWC distribution ----------------------------------------------------------
mu <- c(0, 0)              # mode centered -> bright peak in the middle of the hex
xi <- c(0.3, 0.4, 0.3)    # (xi1, xi2, rho): mild, unequal concentrations + mild rho.
                          # This regime yields a ridge whose connected component
                          # through the mode wraps the whole torus (a gentle wave),
                          # rather than a short central arc.

## Density heatmap over the torus [-pi, pi)^2 --------------------------------
nth <- 300
th <- seq(-pi, pi, length.out = nth)
grid <- expand.grid(theta1 = th, theta2 = th)
grid$dens <- d_bwc(x = as.matrix(grid[, c("theta1", "theta2")]), mu = mu, xi = xi)

## Density ridge (connected component through the mode) ----------------------
ridge <- ridge_bwc(mu = mu, xi = xi, subint_1 = 500, subint_2 = 500)
ridge <- data.frame(theta1 = ridge[, 1], theta2 = ridge[, 2])
ridge <- ridge[order(ridge$theta1), ]
# Break the path where it wraps across the +/- pi boundary, so no spurious
# straight segment is drawn across the panel.
ridge$theta2[c(FALSE, abs(diff(ridge$theta2)) > pi)] <- NA

## Interior subplot: viridis heatmap + red ridge, filling the panel ----------
subplot <- ggplot() +
  geom_raster(data = grid, aes(x = theta1, y = theta2, fill = dens)) +
  geom_path(data = ridge, aes(x = theta1, y = theta2),
            color = "red", linewidth = 1.1) +
  scale_fill_viridis_c() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0))

## Assemble the hex sticker --------------------------------------------------
dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)
sticker(
  subplot = subplot,
  s_x = 1, s_y = 0.92, s_width = 1.45, s_height = 1.45,
  package = "ridgetorus",
  p_x = 1, p_y = 1.52, p_size = 15, p_color = "white",
  h_fill = "#440154", h_color = "#440154",       # viridis dark purple
  url = "github.com/egarpor/ridgetorus", u_size = 3, u_color = "white",
  dpi = 320,
  filename = "man/figures/logo.png"
)

cat("Wrote man/figures/logo.png\n")
