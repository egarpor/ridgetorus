
#
# Hex sticker logo for the ridgetorus package.
#
# The interior is the package's signature object: a bivariate sine von Mises
# (BvM) density on the torus drawn as a viridis heatmap, with its density ridge
# overlaid in red. Low, *unequal* concentrations make it disperse and bend the
# ridge into a curve; a mild dependence keeps it unimodal.
#
# Re-run from the package root with:  Rscript logo/logo.R
# Output:  logo/logo.png (master) and man/figures/logo.png (shipped mirror)
#

library(ridgetorus)
library(ggplot2)
library(hexSticker)
library(magick)

## ---- Shared logo standard (identical across egarpor packages) -------------
# Aller_Rg is bundled with (and auto-registered by) hexSticker, so no
# showtext/font_add setup is needed for the sticker() idiom.
FONT    <- "Aller_Rg"   # typeface for the package name and the GitHub URL
P_SIZE  <- 31.2         # package-name size (shared across packages)
U_SIZE  <- 9.0          # GitHub URL size (large enough to read)
U_X     <- 1.00         # GitHub URL position: along the lower-right hex edge
U_Y     <- 0.08
U_ANGLE <- 30
H_SIZE  <- 1.5          # hexagon border thickness
DPI     <- 600

## Bivariate sine von Mises distribution -------------------------------------
mu    <- c(0, 0)          # mode centered -> bright peak in the middle of the hex
kappa <- c(1, 1.6, 0.7)   # (kappa1, kappa2, lambda): low, *unequal* concentrations
                          # (disperse) with mild dependence. Unimodal because
                          # lambda^2 < kappa1 * kappa2; the unequal concentrations
                          # bend the ridge into a gentle curve.

## Density heatmap over the torus [-pi, pi)^2 --------------------------------
nth <- 300
th <- seq(-pi, pi, length.out = nth)
grid <- expand.grid(theta1 = th, theta2 = th)
grid$dens <- d_bvm(x = as.matrix(grid[, c("theta1", "theta2")]), mu = mu,
                   kappa = kappa)

## Density ridge (connected component through the mode) ----------------------
ridge_raw <- ridge_bvm(mu = mu, kappa = kappa, subint_1 = 500, subint_2 = 500)
ridge <- data.frame(theta1 = ridge_raw[, 1], theta2 = ridge_raw[, 2])
ridge <- ridge[order(ridge$theta1), ]
# Break the path where it wraps across the +/- pi boundary, so no spurious
# straight segment is drawn across the panel.
ridge$theta2[c(FALSE, abs(diff(ridge$theta2)) > pi)] <- NA

## Data points and their orthogonal projections onto the ridge ---------------
# A sample from the model, each point joined to its closest ridge point in the
# torus metric -- the ridge-projection idea shown in show_ridge_pca().
set.seed(15)
data_pts <- r_bvm(n = 20, mu = mu, kappa = kappa)
# Keep the wordmark area (top of the hex) clear: relocate any sampled point that
# would sit under the "ridgetorus" text.
in_name <- function(p) abs(p[1]) < 1.55 & p[2] > 1.05
for (i in which(apply(data_pts, 1, in_name))) {
  repeat {
    q <- r_bvm(n = 1, mu = mu, kappa = kappa)
    if (!in_name(q)) { data_pts[i, ] <- q; break }
  }
}
# Balance the (otherwise all-below) scatter: move three points above the ridge on
# the right, spread along it at distinct feet so their projections do not overlap.
ridge2_at <- function(t1) ridge_raw[which.min(abs(ridge_raw[, 1] - t1)), 2]
balance_t1  <- c(0.62, 1.20, 1.78)   # distinct projection feet along the ridge
balance_off <- c(0.52, 0.46, 0.42)   # heights above the ridge (clear of the name)
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
# Drop any projection segment that would wrap across the +/- pi boundary.
segs <- segs[abs(segs$x - segs$xend) < pi & abs(segs$y - segs$yend) < pi, ]

## Interior subplot: heatmap + red ridge + data points and projections -------
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

## Assemble the hex sticker --------------------------------------------------
# The heatmap is a square, so we scale it up (s_width) until its edges fall
# outside the hexagon and then crop to the hex, so the square outline is never
# visible inside the sticker. hexSticker does not clip an oversized subplot, so
# we render (a) the full-bleed sticker and (b) a hexagon-only silhouette with the
# same geometry, then use (b) as an alpha mask to crop (a) cleanly to the hexagon
# while keeping the wordmark, URL and border.
render <- function(subplot, s_width, pkg, url, file) {
  sticker(
    subplot = subplot, s_x = 1, s_y = 1.0,
    s_width = s_width, s_height = s_width,
    package = pkg, p_x = 1, p_y = 1.52, p_size = P_SIZE,
    p_color = "white", p_family = FONT,
    h_fill = "#440154", h_color = "#5ec962", h_size = H_SIZE,  # viridis purple + green border
    url = url, u_x = U_X, u_y = U_Y, u_angle = U_ANGLE,
    u_size = U_SIZE, u_color = "white", u_family = FONT,
    dpi = DPI, filename = file
  )
  invisible(file)
}

full_f <- tempfile(fileext = ".png")
render(subplot, 2.1, "ridgetorus", "github.com/egarpor/ridgetorus", full_f)

blank  <- image_blank(4, 4, color = "none")     # transparent -> pure hex silhouette
btmp   <- tempfile(fileext = ".png")
image_write(blank, btmp, format = "png")
mask_f <- tempfile(fileext = ".png")
render(btmp, 0.1, "", "", mask_f)               # no wordmark/url on the mask

sticker_img <- image_composite(image_read(full_f), image_read(mask_f),
                               operator = "CopyOpacity")

dir.create("logo", showWarnings = FALSE)
dir.create("man/figures", recursive = TRUE, showWarnings = FALSE)
image_write(sticker_img, "logo/logo.png", format = "png")
image_write(sticker_img, "man/figures/logo.png", format = "png")

cat("Wrote logo/logo.png and man/figures/logo.png\n")
