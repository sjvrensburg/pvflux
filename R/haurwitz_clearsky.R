#' Haurwitz Clear Sky Model
#'
#' @description
#' Calculate clear-sky GHI using the Haurwitz clear sky model.
#'
#' This model estimates clear-sky global horizontal irradiance (GHI) as a function
#' of solar zenith angle only. The Haurwitz model is one of the simpler clear-sky
#' models and requires only solar position data.
#'
#' Based on the implementation in pvlib-python, which follows the Haurwitz clear sky
#' model as presented in Haurwitz (1945, 1946). A report on clear sky models found
#' the Haurwitz model to have the best performance in terms of average monthly error
#' among models which require only zenith angle.
#'
#' The model is valid for solar zenith angles less than 90 degrees.
#'
#' @param time Timestamps as POSIXct, POSIXlt, character, or numeric. If a timezone
#'   is specified, times are internally converted to UTC for solar position
#'   calculations and returned in the original timezone. If no timezone is
#'   specified, UTC is assumed. See \code{\link{time_utils}} for details.
#' @param lat Latitude in degrees
#' @param lon Longitude in degrees
#' @param altitude Altitude above sea level in meters. Default: 0 (not used in Haurwitz model)
#' @param dni_extra Extraterrestrial normal irradiance (W/m^2). If NULL (default),
#'   calculated using Spencer (1971) formula with solar_constant.
#' @param solar_constant Solar constant (W/m^2). Default: 1366.1. Only used if
#'   dni_extra is NULL.
#' @param min_cos_zenith Minimum cosine of zenith angle for calculations.
#'   Default: 0.065 (equivalent to 86.3 degrees)
#'
#' @return Data frame with columns:
#' \describe{
#'   \item{time}{Input timestamps}
#'   \item{ghi_clearsky}{Clear-sky global horizontal irradiance (W/m^2)}
#'   \item{dni_clearsky}{Clear-sky direct normal irradiance (W/m^2)}
#'   \item{dhi_clearsky}{Clear-sky diffuse horizontal irradiance (W/m^2)}
#'   \item{zenith}{Solar zenith angle (degrees)}
#'   \item{airmass}{Relative optical airmass (dimensionless)}
#' }
#'
#' @references
#' Haurwitz, B. (1945). "Insolation in Relation to Cloudiness and Cloud
#' Density," Journal of Meteorology, vol. 2, pp. 154-166.
#'
#' Haurwitz, B. (1946). "Insolation in Relation to Cloud Type," Journal of
#' Meteorology, vol. 3, pp. 123-124.
#'
#' Reno, M. J., Hansen, C. W., and Stein, J. S. (2012). "Global Horizontal
#' Irradiance Clear Sky Models: Implementation and Analysis", Sandia National
#' Laboratories, SAND2012-2389.
#'
#' @examples
#' \dontrun{
#' time <- seq(as.POSIXct("2026-01-15 06:00", tz = "UTC"),
#'             by = "hour", length.out = 12)
#'
#' # Haurwitz clear-sky estimation
#' result <- haurwitz_clearsky(
#'   time = time,
#'   lat = -30.6279,
#'   lon = 24.0054
#' )
#' }
#'
#' @export
haurwitz_clearsky <- function(
  time,
  lat,
  lon,
  altitude = 0,
  dni_extra = NULL,
  solar_constant = 1366.1,
  min_cos_zenith = 0.065
) {
  # Prepare time: convert to UTC for calculations, store original timezone
  time_info <- prepare_time_utc(time)
  time_utc <- time_info$time_utc
  original_tz <- time_info$original_tz

  # Convert time to Julian Day (using UTC time)
  jd <- JD(time_utc)

  # Calculate sun position (timezone = 0 since we're using UTC)
  sv <- sunvector(jd, lat, lon, 0)
  sp <- sunpos(sv)
  theta_z_deg <- sp[, 2]  # Solar zenith angle in degrees

  # Calculate extraterrestrial radiation if not provided
  if (is.null(dni_extra)) {
    doy <- as.numeric(format(time_utc, "%j"))
    dni_extra <- get_extra_radiation_spencer(doy, solar_constant)
  }

  # Calculate airmass (absolute, pressure-corrected)
  cos_zenith <- pmax(cos(theta_z_deg * pi / 180), min_cos_zenith)
  airmass_relative <- kasten_young_airmass(theta_z_deg)
  # Haurwitz model doesn't use altitude, but we still calculate airmass for consistency
  airmass_absolute <- airmass_relative * atm_pressure_altitude_correction(altitude)

  # Haurwitz model implementation
  # The Haurwitz formula: I = 1098 * cos(z) * exp(-0.059 / cos(z))
  # where z is the zenith angle in radians
  ghi <- numeric(length(time_utc))
  cos_zen_gte_0 <- cos_zenith > 0

  if (any(cos_zen_gte_0)) {
    ghi[cos_zen_gte_0] <- (1098.0 * cos_zenith[cos_zen_gte_0] *
                           exp(-0.059 / cos_zenith[cos_zen_gte_0]))
  }

  # For nighttime, set GHI to 0
  ghi[!cos_zen_gte_0] <- 0

  # Calculate DNI using the simplified exponential atmospheric model
  # This is an approximation based on atmospheric transmission
  tau <- 0.6  # Atmospheric transmission coefficient (typical clear value)
  dni <- numeric(length(time_utc))
  dni[cos_zen_gte_0] <- dni_extra[cos_zen_gte_0] * tau^airmass_absolute[cos_zen_gte_0]
  dni[!cos_zen_gte_0] <- 0

  # Calculate DHI as the difference between GHI and DNI*cos(zenith)
  # This follows the same approach as other clear-sky models
  dhi <- ghi - dni * cos_zenith

  # Ensure non-negative values
  ghi <- pmax(ghi, 0)
  dni <- pmax(dni, 0)
  dhi <- pmax(dhi, 0)

  # Restore original timezone for output
  time_out <- restore_time_tz(time_utc, original_tz)

  data.frame(
    time = time_out,
    ghi_clearsky = ghi,
    dni_clearsky = dni,
    dhi_clearsky = dhi,
    zenith = theta_z_deg,
    airmass = airmass_absolute
  )
}