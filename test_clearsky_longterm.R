#!/usr/bin/env Rscript
# Test script for long-term clear-sky pipeline with automatic turbidity lookup

library(pvflux)

cat(strrep("=", 70), "\n")
cat("Test: Long-term clear-sky power simulation with auto turbidity lookup\n")
cat(strrep("=", 70), "\n\n")

# De Aar location
lat <- -30.6279
lon <- 24.0054
altitude <- 1233

# Create 1-year hourly time series
cat("Creating 1-year hourly time series...\n")
time_year <- seq(
  as.POSIXct("2026-01-01 00:00", tz = "Africa/Johannesburg"),
  as.POSIXct("2026-12-31 23:00", tz = "Africa/Johannesburg"),
  by = "hour"
)
cat("  Length:", length(time_year), "hours\n\n")

# Create simple weather data (not realistic, just for testing)
# Temperature: seasonal variation
doy <- as.numeric(format(time_year, "%j"))
T_air <- 20 + 10 * sin((doy - 15) * 2 * pi / 365)  # 10-30°C range

# Wind: constant for simplicity
wind <- rep(2.5, length(time_year))

cat("Running pv_clearsky_dc_pipeline() with automatic turbidity lookup...\n")
cat("  (This will use lookup_linke_turbidity() internally)\n\n")

# Run clear-sky DC pipeline with default linke_turbidity = NULL
system.time({
  result_dc <- pv_clearsky_dc_pipeline(
    time = time_year,
    lat = lat,
    lon = lon,
    T_air = T_air,
    wind = wind,
    tilt = 30,
    azimuth = 0,
    altitude = altitude,
    transposition_model = "haydavies",
    cell_temp_model = "skoplaki",
    P_dc0 = 44880 * 230  # Total DC capacity (W)
  )
})

cat("\nResults summary:\n")
cat("  Rows:", nrow(result_dc), "\n")
cat("  Columns:", ncol(result_dc), "\n")
cat("  Column names:", paste(head(names(result_dc), 10), collapse = ", "), "...\n\n")

# Analyze seasonal variation in clear-sky irradiance
monthly_stats <- aggregate(
  result_dc[, c("ghi_clearsky", "P_dc")],
  by = list(Month = format(result_dc$time, "%B")),
  FUN = function(x) c(mean = mean(x), max = max(x))
)

cat("Monthly average clear-sky GHI and DC power:\n")
month_order <- c("January", "February", "March", "April", "May", "June",
                 "July", "August", "September", "October", "November", "December")
monthly_stats <- monthly_stats[match(month_order, monthly_stats$Month), ]

df_monthly <- data.frame(
  Month = monthly_stats$Month,
  Avg_GHI = round(monthly_stats$ghi_clearsky[, "mean"]),
  Max_GHI = round(monthly_stats$ghi_clearsky[, "max"]),
  Avg_P_dc = round(monthly_stats$P_dc[, "mean"] / 1000, 1),  # Convert to kW
  Max_P_dc = round(monthly_stats$P_dc[, "max"] / 1000, 1)
)
print(df_monthly)

cat("\nAnnual statistics:\n")
cat("  Total clear-sky energy:", round(sum(result_dc$P_dc) / 1e6, 1), "MWh\n")
cat("  Mean daily energy:", round(sum(result_dc$P_dc) / 365 / 1000, 1), "kWh\n")
cat("  Peak power:", round(max(result_dc$P_dc) / 1000, 1), "kW\n")
cat("  Capacity factor:", round(mean(result_dc$P_dc) / max(result_dc$P_dc) * 100, 1), "%\n\n")

cat(strrep("=", 70), "\n")
cat("Test: Compare with fixed turbidity\n")
cat(strrep("=", 70), "\n\n")

# Sample one month to compare database lookup vs fixed turbidity
time_month <- seq(
  as.POSIXct("2026-06-01 00:00", tz = "Africa/Johannesburg"),
  as.POSIXct("2026-06-30 23:00", tz = "Africa/Johannesburg"),
  by = "hour"
)

T_air_month <- rep(20, length(time_month))
wind_month <- rep(2.5, length(time_month))

cat("Comparing June 2026: Database lookup vs. Fixed TL=3.0\n\n")

# Database lookup (default)
result_db <- pv_clearsky_dc_pipeline(
  time = time_month,
  lat = lat,
  lon = lon,
  T_air = T_air_month,
  wind = wind_month,
  tilt = 30,
  azimuth = 0,
  altitude = altitude,
  P_dc0 = 230
)

# Fixed turbidity
result_fixed <- pv_clearsky_dc_pipeline(
  time = time_month,
  lat = lat,
  lon = lon,
  T_air = T_air_month,
  wind = wind_month,
  tilt = 30,
  azimuth = 0,
  altitude = altitude,
  linke_turbidity = 3.0,
  P_dc0 = 230
)

cat("Monthly totals:\n")
cat("  Database lookup:", round(sum(result_db$P_dc) / 1000, 1), "kWh\n")
cat("  Fixed TL=3.0:   ", round(sum(result_fixed$P_dc) / 1000, 1), "kWh\n")
cat("  Difference:     ", round((sum(result_db$P_dc) - sum(result_fixed$P_dc)) / 1000, 1), "kWh\n")
cat("  Percent diff:   ", round((sum(result_db$P_dc) / sum(result_fixed$P_dc) - 1) * 100, 2), "%\n\n")

cat("Peak values:\n")
cat("  Database GHI:", round(max(result_db$ghi_clearsky)), "W/m²\n")
cat("  Fixed GHI:   ", round(max(result_fixed$ghi_clearsky)), "W/m²\n")
cat("  Difference:  ", round(max(result_db$ghi_clearsky) - max(result_fixed$ghi_clearsky)), "W/m²\n\n")

cat(strrep("=", 70), "\n")
cat("Test: AC power pipeline with automatic turbidity\n")
cat(strrep("=", 70), "\n\n")

cat("Running pv_clearsky_power_pipeline() for one week...\n\n")

time_week <- seq(
  as.POSIXct("2026-06-15 00:00", tz = "Africa/Johannesburg"),
  as.POSIXct("2026-06-21 23:00", tz = "Africa/Johannesburg"),
  by = "hour"
)

T_air_week <- rep(25, length(time_week))
wind_week <- rep(2.5, length(time_week))

result_ac <- pv_clearsky_power_pipeline(
  time = time_week,
  lat = lat,
  lon = lon,
  T_air = T_air_week,
  wind = wind_week,
  tilt = 30,
  azimuth = 0,
  altitude = altitude,
  P_dc0 = 44880 * 230,
  n_inverters = 20,
  inverter_kw = 500
)

cat("Weekly summary:\n")
cat("  Total DC energy:", round(sum(result_ac$P_dc) / 1e6, 2), "MWh\n")
cat("  Total AC energy:", round(sum(result_ac$P_ac) / 1e6, 2), "MWh\n")
cat("  Inverter efficiency:", round(sum(result_ac$P_ac) / sum(result_ac$P_dc) * 100, 2), "%\n")
cat("  Any clipping?", any(result_ac$clipped), "\n\n")

cat(strrep("=", 70), "\n")
cat("All long-term tests completed successfully!\n")
cat(strrep("=", 70), "\n")
