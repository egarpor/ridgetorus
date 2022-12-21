
# Required packages
library(ggmap)
library(lubridate)
library(circular)
library(dplyr)
library(tidyverse)

# Set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Plotting the four areas in the map
load("map_santa_barbara.RData")

# Zone A
begin_lat_A <- 34.35
end_lat_A <- 34.44
begin_lon_A <- -120.25
end_lon_A <- -120.1

# Zone B
begin_lat_B <- 34.30
end_lat_B <- 34.4
begin_lon_B <- -119.88
end_lon_B <- -119.73

# Zone C
begin_lat_C <- 34.04
end_lat_C <- 34.12
begin_lon_C <- -120.07
end_lon_C <- -119.92

# Zone D
begin_lat_D <- 33.85
end_lat_D <- 33.95
begin_lon_D <- -119.95
end_lon_D <- -119.82

# List individual RDatas
files <- list.files(path = "../data-raw/", pattern = "*.RData",
                    full.names = TRUE, recursive = FALSE)

# Function that reads over all the files in the directory containing the raw
# data and return the records inside a given area delimited by longitude
# and latitude
extract_data <- function(begin_lat, end_lat, begin_lon, end_lon,
                         begin_y = 2019, end_y = 2022) {

  # Retrieve monthly data
  monthly_data <- lapply(files, function(x) {

    year <- unlist(strsplit(x, "_"))[5]
    if (year >= begin_y & year < end_y) {

      # Load raw data file
      load(x)
      results <- filter(final_dataframe, lon > begin_lon, lon < end_lon,
                        lat > begin_lat, lat < end_lat)

      # Calculate directions and speed
      results$d <- atan2(x = results$u, y = results$v)
      results$speed <- sqrt(results$u^2 + results$v^2)
      return(results)

    }

  })

  # Merge available data
  total_data <- data.frame()
  for (i in seq_along(monthly_data)) {

    total_data <- rbind(total_data, monthly_data[[i]])

  }

  return(total_data)

}

results_A <- extract_data(begin_lat_A, end_lat_A, begin_lon_A, end_lon_A)
coord_A <- unique(na.omit(results_A)[c("lon", "lat")])
rm(results_A)
results_B <- extract_data(begin_lat_B, end_lat_B, begin_lon_B, end_lon_B)
coord_B <- unique(na.omit(results_B)[c("lon", "lat")])
rm(results_B)
results_C <- extract_data(begin_lat_C, end_lat_C, begin_lon_C, end_lon_C)
coord_C <- unique(na.omit(results_C)[c("lon", "lat")])
rm(results_C)
results_D <- extract_data(begin_lat_D, end_lat_D, begin_lon_D, end_lon_D)
coord_D <- unique(na.omit(results_D)[c("lon", "lat")])
rm(results_D)
gc()

png("figures_app/map.png", width = 7, height = 7, units = "in",
    res = 200, bg = "transparent")
ggmap(map.st_barbara) +
  theme(plot.margin = margin(0.1, 0.1, 0, 0, "cm"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  xlab("Longitude") + ylab("Latitude") +
  scale_y_continuous(limits = c(33.7, 34.6)) +
  geom_point(aes(x = lon, y = lat), data = coord_A, colour = 'yellow', size = 1) +
  geom_point(aes(x = lon, y = lat), data = coord_B, colour = 'red', size = 1) +
  geom_point(aes(x = lon, y = lat), data = coord_C, colour = 'pink', size = 1) +
  geom_point(aes(x = lon, y = lat), data = coord_D, colour = 'green', size = 1) +
  annotate(geom = "text", x = -120.18, y = 34.5, label = "A",
                              color = "yellow", size = 7) +
  annotate(geom = "text", x = -119.81, y = 34.45, label = "B",
           color = "red", size = 7) +
  annotate(geom = "text", x = -119.99, y = 34.18, label = "C",
           color = "pink", size = 7) +
  annotate(geom = "text", x = -119.886, y = 33.81, label = "D",
           color = "green", size = 7)
dev.off()

knitr::plot_crop("figures_app/map.png")

# Coordinate limits for each area
print("A")
print(paste(min(coord_A$lon), "x", max(coord_A$lon)))
print(paste(min(coord_A$lat), "x", max(coord_A$lat)))

print("B")
print(paste(min(coord_B$lon), "x", max(coord_B$lon)))
print(paste(min(coord_B$lat), "x", max(coord_B$lat)))

print("C")
print(paste(min(coord_C$lon), "x", max(coord_C$lon)))
print(paste(min(coord_C$lat), "x", max(coord_C$lat)))

print("D")
print(paste(min(coord_D$lon), "x", max(coord_D$lon)))
print(paste(min(coord_D$lat), "x", max(coord_D$lat)))

