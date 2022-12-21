
library(ncdf4)
library(lubridate)
library(circular)
library(dplyr)
library(tidyverse)

# This script shows the necessary steps to obtain the data for the Santa Barbara
# analysis. There are two main steps: (1) downloading the raw data;
# (2) filtering the raw data to obtain data only for zones A, B, C, and D with
# the weighted average daily mean

## Obtaining raw data files for the full west coast

# ID for west coast
data_ID <- "HFR/USWC/2km/hourly/RTV/HFRADAR_US_West_Coast_2km_Resolution_Hourly_RTV_best.ncd"

# Data url
data_url <- paste0("http://hfrnet-tds.ucsd.edu/thredds/dodsC/", data_ID)

# Read data structure
data <- nc_open(data_url)

# Retrieve latitudes, longitudes, and time
lat <- ncvar_get(data, "lat")
lon <- ncvar_get(data, "lon")
tim <- date("2012-01-01 00:00:00 UTC") + hours(ncvar_get(data, "time"))

# Trim download region -- server has a limit of 500Mb request
# Check (lat, lon) coordinates at https://cordc.ucsd.edu/projects/mapping/maps/
# to focus on a given region
loc <- "full_los_angeles_coast"
begin_lat <- 33.7
end_lat <- 34.5
begin_lon <- -121
end_lon <- -119

# Get data function -- takes info from global environment and writes on it
get_data <- function() {

  # Obtain the data indexes associated to the given information
  begin_lat_ind <<- which(lat >= begin_lat)[1]
  begin_lon_ind <<- which(lon >= begin_lon)[1]
  begin_tim_ind <<- which(tim >= begin_tim)[1]
  end_lat_ind <<- which(lat <= end_lat)
  end_lon_ind <<- which(lon <= end_lon)
  end_tim_ind <<- which(tim < end_tim)
  end_lat_ind <<- end_lat_ind[length(end_lat_ind)]
  end_lon_ind <<- end_lon_ind[length(end_lon_ind)]
  end_tim_ind <<- end_tim_ind[length(end_tim_ind)]

  # Download sizes
  l_lat <<- end_lat_ind - begin_lat_ind + 1
  l_lon <<- end_lon_ind - begin_lon_ind + 1
  l_tim <<- end_tim_ind - begin_tim_ind + 1
  stopifnot(l_lat > 0)
  stopifnot(l_lon > 0)
  stopifnot(l_tim > 0)

  # Upper bound on the size of one out of two objects to be downloaded
  cat("Data size:", format(object.size(rnorm(l_lat * l_lon * l_tim)),
                           units = "Mb"))

  # Download (u, v), being:
  # u (m/s) = surface_eastward_sea_water_velocity
  # v (m/s) = surface_northward_sea_water_velocity
  u <<- ncvar_get(data, "u",
                  start = c(begin_lon_ind, begin_lat_ind, begin_tim_ind),
                  count = c(l_lon, l_lat, l_tim))
  v <<- ncvar_get(data, "v",
                  start = c(begin_lon_ind, begin_lat_ind, begin_tim_ind),
                  count = c(l_lon, l_lat, l_tim))

}

# Download data in a monthly loop to avoid surpassing 500Mb limit
for (year in 2017:2022) {

  for (month in 1:12) {

    # Skip from August 2022
    if (year == 2022 && month > 7) {

      break

    }

    # Show progress
    cat(paste0("\n", year, " - ", month, "\n"))

    # Begin and end time taking into account new years
    begin_tim <- date(paste(toString(year), toString(month),
                            "01 00:00:00 UTC", sep = "-"))
    if (month != 12) {

      end_tim <- date(paste(toString(year), toString(month + 1),
                            "01 00:00:00 UTC", sep = "-"))

    } else {

      end_tim <- date(paste(toString(year + 1), toString(1),
                            "01 00:00:00 UTC", sep = "-"))

    }

    # Download data
    get_data()

    # Get data dimensions
    lon_length <- dim(u)[1]
    lat_length <- dim(u)[2]
    time_length <- dim(u)[3]

    # Find the indexes associated to dimensions
    lat_aux <- begin_lat_ind:(begin_lat_ind + lat_length - 1)
    lon_aux <- begin_lon_ind:(begin_lon_ind + lon_length - 1)
    time_aux <- begin_tim_ind:(begin_tim_ind + time_length - 1)

    # Join all the cases respecting the original order of u[lon, lat, time]
    join <- merge(x = lon[lon_aux], y = lat[lat_aux])
    final_dataframe <- merge(x = join, y = tim[time_aux], by = NULL)
    colnames(final_dataframe) <- c("lon", "lat", "time")

    # Add the velocities and save the data
    final_dataframe$u <- c(u)
    final_dataframe$v <- c(v)

    # Order data: time, lon, lat
    ord <- order(final_dataframe$time, final_dataframe$lon, final_dataframe$lat)
    final_dataframe <- final_dataframe[ord, ]

    # Save the data
    save(final_dataframe, file = paste(loc, "_", toString(year), "_",
                                       toString(month), ".RData", sep = ""))

  }

}

## Obtaining daily data for zones A, B, C, and D

# List individual RDatas
files <- list.files(pattern = "*.RData", full.names = TRUE, recursive = FALSE)

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
      results <- filter(final_dataframe,
                        lon > begin_lon, lon < end_lon,
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

# Function that calculates the speed-weighted average of the directions over
# periods of certain number of hours
extract_theta <- function(results, hours) {

  # Assign group according to the desired period
  results$group <- cut(results$time, breaks = paste(hours, "hours"))
  results <- results[complete.cases(results), ]
  # Apply the speed-weighted average over the non-NA observations inside the
  # t hour period

  speedw_mean <- function(x) {

    dir_speeds <- x[, c("d", "speed")]
    dir_speeds <- dir_speeds[complete.cases(dir_speeds), ]
    weights <- dir_speeds$speed / sum(dir_speeds$speed)
    circular:::WeightedMeanCircularRad(w = weights, x = dir_speeds$d)

  }

  theta <- results[, c("group", "d", "speed")] %>% 
    group_by(group) %>%
    do(data.frame(val = speedw_mean(.)))

  return(theta)

}

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

# Obtain the results for every area and their corresponding daily directions
results_A <- extract_data(begin_lat_A, end_lat_A, begin_lon_A, end_lon_A)
A <- extract_theta(results_A, hours = 24)

results_B <- extract_data(begin_lat_B, end_lat_B, begin_lon_B, end_lon_B)
B <- extract_theta(results_B, hours = 24)

results_C <- extract_data(begin_lat_C, end_lat_C, begin_lon_C, end_lon_C)
C <- extract_theta(results_C, hours = 24)

results_D <- extract_data(begin_lat_D, end_lat_D, begin_lon_D, end_lon_D)
D <- extract_theta(results_D, hours = 24)

# Merge results
df_list <- list(A, B, C, D)
santabarbara <- as.data.frame(na.omit(df_list %>%
                                        reduce(full_join, by = 'group'))[, 2:5])

# Change column names
colnames(santabarbara) <- c('A', 'B', 'C', 'D')

# Save the object
save(list = "santabarbara", file = "santabarbara.rda", compress = "bzip2")

# Cleanup
rm(results_A, results_B, results_C, results_D)
rm(A, B, C, D)
gc()
