
# Hex sticker for the 'ridgetorus' package: a bivariate sine von Mises density
# on the torus as a viridis heatmap with its density ridge overlaid in red, and
# a small sample projected onto the ridge. Built with the package's own d_bvm /
# ridge_bvm / r_bvm, then cropped to the hexagon with an alpha mask.

library(ridgetorus)
library(ggplot2)
library(hexSticker)
library(magick)

# Shared logo standards
font <- "Aller_Rg"
name_size <- 31.2
url_size <- 9.0
url_x <- 1.00
url_y <- 0.08
url_angle <- 30
hex_border <- 1.5
dpi <- 600

# Bivariate sine von Mises density: low, unequal concentrations keep it unimodal
# (lambda^2 < kappa1 * kappa2) but bend the ridge into a gentle curve.
mu <- c(0, 0)
kappa <- c(1, 1.6, 0.7) # (kappa1, kappa2, lambda)

# Density heatmap over the torus
nth <- 300
th <- seq(-pi, pi, length.out = nth)
grid <- expand.grid(theta1 = th, theta2 = th)
grid$dens <- d_bvm(x = as.matrix(grid[, c("theta1", "theta2")]), mu = mu,
                   kappa = kappa)

# Density ridge (connected component through the mode)
ridge_raw <- ridge_bvm(mu = mu, kappa = kappa, subint_1 = 500, subint_2 = 500)
ridge <- data.frame(theta1 = ridge_raw[, 1], theta2 = ridge_raw[, 2])
ridge <- ridge[order(ridge$theta1), ]
# Break the path where it wraps across the +/- pi boundary.
ridge$theta2[c(FALSE, abs(diff(ridge$theta2)) > pi)] <- NA

# Sample points and their orthogonal projections onto the ridge
# A model sample, each point joined to its nearest ridge point (torus metric).
set.seed(15)
data_pts <- r_bvm(n = 20, mu = mu, kappa = kappa)
# Keep the wordmark area clear: relocate any point that sits under the name.
in_name <- function(p) abs(p[1]) < 1.55 & p[2] > 1.05
for (i in which(apply(data_pts, 1, in_name))) {
  repeat {
    q <- r_bvm(n = 1, mu = mu, kappa = kappa)
    if (!in_name(q)) { data_pts[i, ] <- q; break }
  }
}
# Move three points above the ridge (right side) so projections don't overlap.
ridge2_at <- function(t1) ridge_raw[which.min(abs(ridge_raw[, 1] - t1)), 2]
balance_t1 <- c(0.62, 1.20, 1.78) # projection feet along the ridge
balance_off <- c(0.52, 0.46, 0.42) # heights above the ridge
for (k in 1:3) {
  i <- c(8, 11, 18)[k]
  data_pts[i, 1] <- balance_t1[k]
  data_pts[i, 2] <- ridge2_at(balance_t1[k]) + balance_off[k]
}
wrap <- function(a) ((a + pi) %% (2 * pi)) - pi
proj_idx <- apply(data_pts, 1, function(th)
  which.min(wrap(ridge_raw[, 1] - th[1])^2 + wrap(ridge_raw[, 2] - th[2])^2))
data_pts <- data.frame(theta1 = data_pts[, 1], theta2 = data_pts[, 2])
proj_pts <- data.frame(theta1 = ridge_raw[proj_idx, 1],
                       theta2 = ridge_raw[proj_idx, 2])
segs <- data.frame(x = data_pts$theta1, y = data_pts$theta2,
                   xend = proj_pts$theta1, yend = proj_pts$theta2)
# Drop any projection segment that wraps across the +/- pi boundary.
segs <- segs[abs(segs$x - segs$xend) < pi & abs(segs$y - segs$yend) < pi, ]

# Interior subplot: heatmap + red ridge + points and projections
subplot <- ggplot() +
  geom_raster(data = grid, aes(x = theta1, y = theta2, fill = dens)) +
  geom_path(data = ridge, aes(x = theta1, y = theta2),
            color = "red", linewidth = 1.1) +
  geom_segment(data = segs, aes(x = x, y = y, xend = xend, yend = yend),
               color = "white", linewidth = 0.35, alpha = 0.6) +
  geom_point(data = proj_pts, aes(x = theta1, y = theta2),
             color = "white", size = 0.7) +
  geom_point(data = data_pts, aes(x = theta1, y = theta2), shape = 21,
             fill = "white", color = "grey15", size = 1.7, stroke = 0.3) +
  scale_fill_viridis_c() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0))

# Assemble the hex sticker
# The square heatmap is scaled past the hex edges, then cropped with an alpha
# mask: (a) the full sticker and (b) a hex-only silhouette give (a)'s alpha, as
# hexSticker does not clip an oversized subplot.
render <- function(subplot, s_width, pkg, url, file) {
  sticker(
    subplot = subplot, s_x = 1, s_y = 1.0,
    s_width = s_width, s_height = s_width,
    package = pkg, p_x = 1, p_y = 1.52, p_size = name_size,
    p_color = "white", p_family = font,
    h_fill = "#440154", h_color = "#5ec962", h_size = hex_border,
    url = url,
    u_x = url_x, u_y = url_y, u_angle = url_angle, u_size = url_size,
    u_color = "white", u_family = font,
    dpi = dpi, filename = file
  )
  invisible(file)
}

# Render the full-bleed sticker
full_f <- tempfile(fileext = ".png")
render(subplot, 2.1, "ridgetorus", "github.com/egarpor/ridgetorus", full_f)

# Render a hexagon-only silhouette to use as the alpha mask
blank <- image_blank(4, 4, color = "none")
btmp <- tempfile(fileext = ".png")
image_write(blank, btmp, format = "png")
mask_f <- tempfile(fileext = ".png")
render(btmp, 0.1, "", "", mask_f)

# Crop to the hexagon using the silhouette as the alpha channel
sticker_img <- image_composite(image_read(full_f), image_read(mask_f),
                               operator = "CopyOpacity")

# Write the cropped sticker to logo/ and man/figures/
dir.create("logo", showWarnings = FALSE)
dir.create("man/figures", recursive = TRUE, showWarnings = FALSE)
image_write(sticker_img, "logo/logo.png", format = "png")
image_write(sticker_img, "man/figures/logo.png", format = "png")
