#' Hay and Davies Anisotropic Transposition Model
#'
#' @description
#' Convert horizontal irradiance to plane-of-array (POA) irradiance using the
#' Hay and Davies (1980) anisotropic sky model.
#'
#' The Hay and Davies model determines the total POA irradiance by combining:
#' \itemize{
#'   \item Beam component: Direct normal irradiance projected onto the tilted surface
#'   \item Sky diffuse component: Anisotropic model with circumsolar and isotropic parts
#'   \item Ground diffuse component: Isotropic reflection from the ground
#' }
#'
#' The sky diffuse irradiance is calculated as:
#' \deqn{I_{d} = DHI \left( A\cdot R_b + (1 - A) \left(\frac{1 + \cos\beta}{2}\right ) \right )}
#'
#' where \eqn{A = DNI / I_{extra}} is the anisotropy index, \eqn{R_b} is the
#' projection ratio (cosine of angle of incidence to cosine of zenith angle),
#' and \eqn{\beta} is the tilt angle of the array.
#'
#' @param time Timestamps as POSIXct, POSIXlt, character, or numeric. If a timezone
#'   is specified, times are internally converted to UTC for solar position
#'   calculations and returned in the original timezone. If no timezone is
#'   specified, UTC is assumed. See \code{\link{time_utils}} for details.
#' @param lat Latitude in degrees
#' @param lon Longitude in degrees
#' @param GHI Global horizontal irradiance (W/m^2)
#' @param DNI Direct normal irradiance (W/m^2)
#' @param DHI Diffuse horizontal irradiance (W/m^2)
#' @param tilt Panel tilt angle (degrees from horizontal)
#' @param azimuth Panel azimuth (degrees, 0 = north)
#' @param albedo Ground albedo (default 0.2)
#' @param min_cos_zenith Minimum value of cos(zenith) for Rb calculation. Default: 0.01745
#' @param solar_constant The solar constant (W/m^2). Default: 1366.1
#'
#' @return Data frame with columns:
#' \describe{
#'   \item{time}{Input timestamps}
#'   \item{GHI}{Input global horizontal irradiance (W/m^2)}
#'   \item{DNI}{Input direct normal irradiance (W/m^2)}
#'   \item{DHI}{Input diffuse horizontal irradiance (W/m^2)}
#'   \item{poa_global}{Total plane-of-array irradiance (W/m^2)}
#'   \item{poa_beam}{Beam (direct) component on tilted surface (W/m^2)}
#'   \item{poa_sky_diffuse}{Sky diffuse component on tilted surface (W/m^2)}
#'   \item{poa_ground_diffuse}{Ground reflected component on tilted surface (W/m^2)}
#'   \item{poa_diffuse}{Total diffuse (sky + ground) on tilted surface (W/m^2)}
#'   \item{zenith}{Solar zenith angle (degrees)}
#'   \item{azimuth}{Solar azimuth angle (degrees)}
#'   \item{incidence}{Angle of incidence on panel (degrees)}
#'   \item{ai}{Anisotropy index}
#'   \item{rb}{Projection ratio}
#' }
#'
#' @references
#' Hay, J. E., & Davies, J. A. (1980). Calculations of the solar radiation
#' incident on an inclined surface. In J. E. Hay & T. K. Won (Eds.),
#' Proc. of First Canadian Solar Radiation Data Workshop (pp. 59).
#' Ministry of Supply and Services, Canada.
#'
#' Loutzenhiser, P. G., et al. (2007). Empirical validation of models to
#' compute solar irradiance on inclined surfaces for building energy simulation.
#' Solar Energy, 81, 254-267. \doi{10.1016/j.solener.2006.03.009}
#'
#' @examples
#' \dontrun{
#' time <- seq(as.POSIXct("2026-01-15 06:00", tz = "UTC"),
#'             by = "hour", length.out = 12)
#' GHI <- c(50, 200, 450, 700, 850, 950, 1000, 950, 850, 700, 450, 200)
#'
#' # First decompose GHI to get DNI and DHI
#' decomposed <- erbs_decomposition(time, lat = -30.6279, lon = 24.0054, GHI = GHI)
#'
#' # Then apply Hay-Davies transposition
#' result <- haydavies_transposition(
#'   time = time,
#'   lat = -30.6279,
#'   lon = 24.0054,
#'   GHI = GHI,
#'   DNI = decomposed$DNI,
#'   DHI = decomposed$DHI,
#'   tilt = 20,
#'   azimuth = 0
#' )
#' }
#'
#' @export
#'
haydavies_transposition <- function(
  time,
  lat,
  lon,
  GHI,
  DNI,
  DHI,
  tilt,
  azimuth,
  albedo = 0.2,
  min_cos_zenith = 0.01745,
  solar_constant = 1366.1
) {
  # Prepare time: convert to UTC for calculations, store original timezone
  time_info <- prepare_time_utc(time)
  time_utc <- time_info$time_utc
  original_tz <- time_info$original_tz

  stopifnot(length(time_utc) == length(GHI))
  stopifnot(length(time_utc) == length(DNI))
  stopifnot(length(time_utc) == length(DHI))

  # Convert time to Julian Day (using UTC time)
  jd <- JD(time_utc)

  # Calculate sun vector and position (timezone = 0 since we're using UTC)
  sv <- sunvector(jd, lat, lon, 0)
  sp <- sunpos(sv)

  # Extract zenith and azimuth
  theta_z_deg <- sp[, 2]  # Solar zenith angle in degrees
  sun_az_deg <- sp[, 1]   # Solar azimuth angle in degrees

  # =========================================================================
  # Calculate Angle of Incidence (AOI)
  # =========================================================================

  # Panel tilt and azimuth in radians
  beta <- tilt * pi / 180
  panel_az <- azimuth * pi / 180

  # Solar zenith and azimuth in radians
  theta_z <- theta_z_deg * pi / 180
  sun_az <- sun_az_deg * pi / 180

  # Calculate AOI projection (cosine of angle of incidence)
  # This is the dot product of panel normal and solar position unit vectors
  cos_aoi <- (cos(beta) * cos(theta_z) +
                sin(beta) * sin(theta_z) * cos(sun_az - panel_az))
  cos_aoi <- pmax(-1, pmin(1, cos_aoi))  # Clip to valid range

  # Angle of incidence in degrees
  aoi_deg <- acos(cos_aoi) * 180 / pi

  # =========================================================================
  # Calculate Projection Ratio (Rb)
  # =========================================================================

  cos_z <- cos(theta_z)
  cos_aoi <- pmax(0, cos_aoi)  # Negative values when sun is behind panel
  rb <- cos_aoi / pmax(cos_z, min_cos_zenith)

  # =========================================================================
  # Calculate Extraterrestrial Radiation and Anisotropy Index
  # =========================================================================

  doy <- as.numeric(format(time_utc, "%j"))
  dni_extra <- get_extra_radiation_spencer(doy, solar_constant)

  # Anisotropy index: ratio of beam irradiance to extraterrestrial irradiance
  ai <- DNI / dni_extra
  ai <- pmax(0, pmin(1, ai))  # Constrain to [0, 1]

  # =========================================================================
  # Hay and Davies (1980) Sky Diffuse Model
  # =========================================================================

  # Isotropic component term
  term1 <- 1 - ai
  term2 <- 0.5 * (1 + cos(beta))

  # Sky diffuse components
  poa_isotropic <- pmax(0, DHI * term1 * term2)
  poa_circumsolar <- pmax(0, DHI * ai * rb)

  # Total sky diffuse
  poa_sky_diffuse <- poa_isotropic + poa_circumsolar

  # =========================================================================
  # Beam Component (Direct Normal Irradiance projected onto tilted surface)
  # =========================================================================

  # Only beam irradiance when sun is in front of panel
  poa_beam <- DNI * pmax(0, cos_aoi)

  # =========================================================================
  # Ground Reflected Diffuse Component (Isotropic)
  # =========================================================================

  poa_ground_diffuse <- GHI * albedo * 0.5 * (1 - cos(beta))

  # =========================================================================
  # Total POA Irradiance
  # =========================================================================

  poa_diffuse <- poa_sky_diffuse + poa_ground_diffuse
  poa_global <- poa_beam + poa_diffuse

  # Set nighttime values to zero
  nighttime <- cos_z <= 0
  poa_global[nighttime] <- 0
  poa_beam[nighttime] <- 0
  poa_sky_diffuse[nighttime] <- 0
  poa_ground_diffuse[nighttime] <- 0
  poa_diffuse[nighttime] <- 0

  # Restore original timezone for output
  time_out <- restore_time_tz(time_utc, original_tz)

  data.frame(
    time = time_out,
    GHI = GHI,
    DNI = DNI,
    DHI = DHI,
    poa_global = poa_global,
    poa_beam = poa_beam,
    poa_sky_diffuse = poa_sky_diffuse,
    poa_ground_diffuse = poa_ground_diffuse,
    poa_diffuse = poa_diffuse,
    zenith = theta_z_deg,
    azimuth = sun_az_deg,
    incidence = aoi_deg,
    ai = ai,
    rb = rb
  )
}
