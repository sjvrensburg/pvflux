#' Erbs Decomposition of GHI into DNI and DHI
#'
#' @description
#' Estimate DNI and DHI from GHI using the Erbs model (Erbs et al., 1982).
#'
#' The Erbs model estimates the diffuse fraction (DF) from global horizontal
#' irradiance through an empirical relationship between DF and the ratio of
#' GHI to extraterrestrial irradiance (Kt).
#'
#' \deqn{DHI = DF \times GHI}
#' \deqn{DNI = (GHI - DHI) / \cos(\theta_z)}
#'
#' where \eqn{\theta_z} is the solar zenith angle.
#'
#' @param time Timestamps as POSIXct, POSIXlt, character, or numeric. If a timezone
#'   is specified, times are internally converted to UTC for solar position
#'   calculations and returned in the original timezone. If no timezone is
#'   specified, UTC is assumed. See \code{\link{time_utils}} for details.
#' @param lat Latitude in degrees
#' @param lon Longitude in degrees
#' @param GHI Global horizontal irradiance (W/m^2)
#' @param min_cos_zenith Minimum value of cos(zenith) to allow when calculating
#'   global clearness index kt. Equivalent to zenith = 86.273 degrees.
#'   Default: 0.065
#' @param max_zenith Maximum value of zenith to allow in DNI calculation.
#'   DNI will be set to 0 for times with zenith values greater than max_zenith.
#'   Default: 87
#' @param solar_constant The solar constant (W/m^2). Default: 1366.1
#'
#' @return Data frame with columns:
#' \describe{
#'   \item{time}{Input timestamps}
#'   \item{GHI}{Input global horizontal irradiance (W/m^2)}
#'   \item{DNI}{Direct normal irradiance (W/m^2)}
#'   \item{DHI}{Diffuse horizontal irradiance (W/m^2)}
#'   \item{kt}{Clearness index [unitless]}
#'   \item{zenith}{Solar zenith angle (degrees)}
#' }
#'
#' @references
#' Erbs, D. G., Klein, S. A., & Duffie, J. A. (1982). Estimation of the
#' diffuse radiation fraction for hourly, daily and monthly-average global
#' radiation. Solar Energy, 28(4), 293-302.
#' \doi{10.1016/0038-092X(82)90302-4}
#'
#' @examples
#' \dontrun{
#' time <- seq(as.POSIXct("2026-01-15 06:00", tz = "UTC"),
#'             by = "hour", length.out = 12)
#' GHI <- c(50, 200, 450, 700, 850, 950, 1000, 950, 850, 700, 450, 200)
#'
#' result <- erbs_decomposition(
#'   time = time,
#'   lat = -30.6279,
#'   lon = 24.0054,
#'   GHI = GHI
#' )
#' }
#'
#' @export
#'
erbs_decomposition <- function(
  time,
  lat,
  lon,
  GHI,
  min_cos_zenith = 0.065,
  max_zenith = 87,
  solar_constant = 1366.1
) {
  # Prepare time: convert to UTC for calculations, store original timezone
  time_info <- prepare_time_utc(time)
  time_utc <- time_info$time_utc
  original_tz <- time_info$original_tz

  stopifnot(length(time_utc) == length(GHI))

  # Convert time to Julian Day (using UTC time)
  jd <- JD(time_utc)

  # Calculate sun position (timezone = 0 since we're using UTC)
  sv <- sunvector(jd, lat, lon, 0)
  sp <- sunpos(sv)
  theta_z_deg <- sp[, 2]  # Solar zenith angle in degrees

  # Get day of year (from UTC time)
  doy <- as.numeric(format(time_utc, "%j"))

  # Calculate extraterrestrial radiation (normal to sun)
  dni_extra <- get_extra_radiation_spencer(doy, solar_constant)

  # Calculate clearness index
  kt <- clearness_index(GHI, theta_z_deg, dni_extra, min_cos_zenith)

  # Calculate diffuse fraction using Erbs model
  # For Kt <= 0.22
  df <- 1 - 0.09 * kt

  # For 0.22 < Kt <= 0.8
  df <- ifelse(kt > 0.22 & kt <= 0.8,
               0.9511 - 0.1604 * kt + 4.388 * kt^2 -
                 16.638 * kt^3 + 12.336 * kt^4,
               df)

  # For Kt > 0.8
  df <- ifelse(kt > 0.8, 0.165, df)

  # Calculate DHI
  dhi <- df * GHI

  # Calculate DNI
  cos_zenith <- cos(theta_z_deg * pi / 180)
  dni <- (GHI - dhi) / cos_zenith

  # Handle bad values
  bad_values <- (theta_z_deg > max_zenith) | (GHI < 0) | (dni < 0)
  dni <- ifelse(bad_values, 0, dni)
  dhi <- ifelse(bad_values, GHI, dhi)

  # Restore original timezone for output
  time_out <- restore_time_tz(time_utc, original_tz)

  data.frame(
    time = time_out,
    GHI = GHI,
    DNI = dni,
    DHI = dhi,
    kt = kt,
    zenith = theta_z_deg
  )
}

#' Extraterrestrial Radiation (Spencer Method)
#'
#' @description
#' Determine extraterrestrial radiation from day of year using the
#' Spencer (1971) method.
#'
#' @param doy Day of year (numeric or vector)
#' @param solar_constant The solar constant (W/m^2). Default: 1366.1
#'
#' @return Extraterrestrial radiation normal to the sun (W/m^2)
#'
#' @references
#' Spencer, J. W. (1971). Fourier series representation of the sun. Search, 2, 172.
#'
#' @keywords internal
#'
get_extra_radiation_spencer <- function(doy, solar_constant = 1366.1) {
  # Day angle (Spencer method uses offset=1)
  B <- (2 * pi / 365) * (doy - 1)

  # Earth-Sun distance correction factor (Spencer, 1971)
  RoverR0sqrd <- (1.00011 + 0.034221 * cos(B) + 0.00128 * sin(B) +
                   0.000719 * cos(2 * B) + 0.000077 * sin(2 * B))

  # Extraterrestrial radiation
  Ea <- solar_constant * RoverR0sqrd

  return(Ea)
}

#' Clearness Index
#'
#' @description
#' Calculate the clearness index, the ratio of global to extraterrestrial
#' irradiance on a horizontal plane.
#'
#' @param ghi Global horizontal irradiance (W/m^2)
#' @param solar_zenith Solar zenith angle in degrees
#' @param extra_radiation Extraterrestrial radiation (W/m^2)
#' @param min_cos_zenith Minimum value of cos(zenith) to allow. Default: 0.065
#' @param max_clearness_index Maximum value of the clearness index. Default: 2.0
#'
#' @return Clearness index [unitless]
#'
#' @references
#' Maxwell, E. L. (1987). A Quasi-Physical Model for Converting Hourly
#' Global Horizontal to Direct Normal Insolation. Technical Report No.
#' SERI/TR-215-3087, Solar Energy Research Institute.
#'
#' @keywords internal
#'
clearness_index <- function(ghi, solar_zenith, extra_radiation,
                             min_cos_zenith = 0.065,
                             max_clearness_index = 2.0) {
  cos_zenith <- cos(solar_zenith * pi / 180)
  I0h <- extra_radiation * pmax(cos_zenith, min_cos_zenith)

  kt <- ghi / I0h
  kt <- pmax(kt, 0)
  kt <- pmin(kt, max_clearness_index)

  return(kt)
}
