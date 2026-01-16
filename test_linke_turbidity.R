#!/usr/bin/env Rscript
# Test script for Linke turbidity lookup functionality

library(pvflux)

cat(strrep("=", 70), "\n")
cat("Test 1: Verify LinkeTurbidities.h5 file exists\n")
cat(strrep("=", 70), "\n\n")

db_path <- linke_turbidity_filepath()
cat("Database path:", db_path, "\n")
cat("File exists:", file.exists(db_path), "\n")
cat("File size:", file.info(db_path)$size / 1024^2, "MB\n\n")

cat(strrep("=", 70), "\n")
cat("Test 2: Lookup turbidity for De Aar, monthly values\n")
cat(strrep("=", 70), "\n\n")

# De Aar location
lat <- -30.6279
lon <- 24.0054

# Monthly time series
time_monthly <- seq(
  as.POSIXct("2026-01-15", tz = "UTC"),
  by = "month",
  length.out = 12
)

# Lookup turbidity with interpolation
tl_interp <- lookup_linke_turbidity(
  time = time_monthly,
  lat = lat,
  lon = lon,
  interp_turbidity = TRUE
)

# Lookup turbidity without interpolation
tl_monthly <- lookup_linke_turbidity(
  time = time_monthly,
  lat = lat,
  lon = lon,
  interp_turbidity = FALSE
)

df_monthly <- data.frame(
  Month = format(time_monthly, "%B"),
  TL_Interp = round(tl_interp, 3),
  TL_Monthly = round(tl_monthly, 3)
)

print(df_monthly)

cat("\nMean turbidity (interpolated):", round(mean(tl_interp), 3), "\n")
cat("Range:", round(min(tl_interp), 3), "-", round(max(tl_interp), 3), "\n\n")

cat(strrep("=", 70), "\n")
cat("Test 3: Lookup turbidity for multiple locations\n")
cat(strrep("=", 70), "\n\n")

# Test multiple locations at same time
time_single <- as.POSIXct("2026-06-15", tz = "UTC")
locations <- data.frame(
  Name = c("De Aar (ZA)", "New York (US)", "London (UK)", "Beijing (CN)"),
  Lat = c(-30.6279, 40.7128, 51.5074, 39.9042),
  Lon = c(24.0054, -74.0060, -0.1278, 116.4074)
)

tl_locations <- lookup_linke_turbidity(
  time = rep(time_single, nrow(locations)),
  lat = locations$Lat,
  lon = locations$Lon
)

locations$Turbidity <- round(tl_locations, 3)
print(locations)
cat("\n")

cat(strrep("=", 70), "\n")
cat("Test 4: Compare with simple_linke_turbidity()\n")
cat(strrep("=", 70), "\n\n")

# Compare database lookup with simple estimation
tl_database <- lookup_linke_turbidity(time_monthly, lat, lon)
tl_simple <- simple_linke_turbidity(
  time_monthly,
  location_type = "rural",
  hemisphere = "south"
)

df_compare <- data.frame(
  Month = format(time_monthly, "%B"),
  Database = round(tl_database, 3),
  Simple = round(tl_simple, 3),
  Diff = round(tl_database - tl_simple, 3)
)

print(df_compare)
cat("\nMean absolute difference:", round(mean(abs(df_compare$Diff)), 3), "\n\n")

cat(strrep("=", 70), "\n")
cat("Test 5: Use in ineichen_clearsky() function\n")
cat(strrep("=", 70), "\n\n")

# Create hourly time series for one day
time_hourly <- seq(
  as.POSIXct("2026-06-15 06:00", tz = "Africa/Johannesburg"),
  as.POSIXct("2026-06-15 18:00", tz = "Africa/Johannesburg"),
  by = "hour"
)

# Clear-sky with fixed turbidity
cs_fixed <- ineichen_clearsky(
  time = time_hourly,
  lat = lat,
  lon = lon,
  linke_turbidity = 3.0,
  altitude = 1233
)

# Clear-sky with database lookup
cs_lookup <- ineichen_clearsky(
  time = time_hourly,
  lat = lat,
  lon = lon,
  linke_turbidity = NULL,  # Auto-lookup
  altitude = 1233
)

df_clearsky <- data.frame(
  Time = format(time_hourly, "%H:%M"),
  GHI_Fixed = round(cs_fixed$ghi_clearsky),
  GHI_Lookup = round(cs_lookup$ghi_clearsky),
  Diff = round(cs_lookup$ghi_clearsky - cs_fixed$ghi_clearsky)
)

print(head(df_clearsky, 8))
cat("\nMax GHI (fixed TL=3.0):", round(max(cs_fixed$ghi_clearsky)), "W/m²\n")
cat("Max GHI (database lookup):", round(max(cs_lookup$ghi_clearsky)), "W/m²\n")
cat("Mean difference:", round(mean(abs(df_clearsky$Diff))), "W/m²\n\n")

cat(strrep("=", 70), "\n")
cat("Test 6: Daily variation in turbidity\n")
cat(strrep("=", 70), "\n\n")

# Check turbidity variation throughout one day
time_daily <- seq(
  as.POSIXct("2026-06-15 00:00", tz = "UTC"),
  as.POSIXct("2026-06-15 23:00", tz = "UTC"),
  by = "hour"
)

tl_daily <- lookup_linke_turbidity(time_daily, lat, lon, interp_turbidity = TRUE)

cat("Hourly turbidity for 2026-06-15 at De Aar:\n")
cat("  Min:", round(min(tl_daily), 4), "\n")
cat("  Max:", round(max(tl_daily), 4), "\n")
cat("  Range:", round(max(tl_daily) - min(tl_daily), 4), "\n")
cat("  (Should be nearly constant for daily interpolation)\n\n")

cat(strrep("=", 70), "\n")
cat("All tests completed successfully!\n")
cat(strrep("=", 70), "\n")
