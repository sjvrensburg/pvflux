#' Boland-Ridley Decomposition of GHI into DNI and DHI
#'
#' @description
#' Estimate DNI and DHI from GHI using the Boland-Ridley clearness index model
#' (Boland et al., 2001; Boland & Ridley, 2008).
#'
#' The Boland model estimates the diffuse fraction (DF) from global horizontal
#' irradiance through an empirical relationship between DF and the clearness
#' index (Kt), the ratio of GHI to extraterrestrial irradiance.
#'
#' The Boland model uses a logistic function to fit the entire range of
#' clearness index, unlike other decomposition models that use piecewise
#' polynomial functions:
#' \deqn{DF = \\frac{1}{1 + \\exp\\left(a \\left(k_t - b\\right)\\right)}}
#'
#' @param time POSIXct vector of times (UTC recommended)
#' @param lat Latitude in degrees
#' @param lon Longitude in degrees
#' @param GHI Global horizontal irradiance (W/m^2)
#' @param a_coeff Logistic curve fit coefficient. Default 7.997.
#' @param b_coeff Logistic curve fit coefficient. Default 0.586.
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
#'   \item{df}{Diffuse fraction [unitless]}
#' }
#'
#' @references
#' Boland, J., Ridley, B. (2008). Models of Diffuse Solar Fraction. In:
#' Badescu V. (eds) Modeling Solar Radiation at the Earth's Surface.
#' Springer, Berlin, Heidelberg. \doi{10.1007/978-3-540-77455-6_8}
#'
#' Boland, J., Scott, L., & Luther, M. (2001). Modelling the diffuse fraction
#' of global solar radiation on a horizontal surface. Environmetrics, 12(2),
#' 103-116. \doi{10.1002/1099-095X(200103)12:2<103::AID-ENV447>3.0.CO;2-2}
#'
#' @note
#' Boland diffuse fraction differs from other decomposition algorithms by use
#' of a logistic function to fit the entire range of clearness index, Kt.
#' Parameters `a_coeff` and `b_coeff` are reported in Boland et al. (2001)
#' for different time intervals:
#' \itemize{
#'   \item 15-minute: `a = 8.645` and `b = 0.613`
#'   \item 1-hour: `a = 7.997` and `b = 0.586` (default)
#' }
#'
#' @seealso
#' \code{\\link{erbs_decomposition}} for Erbs decomposition model
#'
#' @examples
#' \dontrun{
#' time <- seq(as.POSIXct("2026-01-15 06:00", tz = "UTC"),
#'             by = "hour", length.out = 12)
#' GHI <- c(50, 200, 450, 700, 850, 950, 1000, 950, 850, 700, 450, 200)
#'
#' result <- boland_decomposition(
#'   time = time,
#'   lat = -30.6279,
#'   lon = 24.0054,
#'   GHI = GHI
#' )
#' }
#'
#' @export
#'
boland_decomposition <- function(
  time,
  lat,
  lon,
  GHI,
  a_coeff = 7.997,
  b_coeff = 0.586,
  min_cos_zenith = 0.065,
  max_zenith = 87,
  solar_constant = 1366.1
) {
  stopifnot(length(time) == length(GHI))

  # Convert time to Julian Day
  jd <- JD(time)

  # Extract timezone from POSIXct object (in hours)
  tz_offset <- as.numeric(format(time[1], "%z")) / 100

  # Calculate sun position
  sv <- sunvector(jd, lat, lon, tz_offset)
  sp <- sunpos(sv)
  theta_z_deg <- sp[, 2]  # Solar zenith angle in degrees

  # Get day of year
  doy <- as.numeric(format(time, "%j"))

  # Calculate extraterrestrial radiation (normal to sun)
  dni_extra <- get_extra_radiation_spencer(doy, solar_constant)

  # Calculate clearness index
  kt <- clearness_index(GHI, theta_z_deg, dni_extra, min_cos_zenith)

  # Boland equation: logistic function for diffuse fraction
  # DF = 1 / (1 + exp(a * (kt - b)))
  df <- 1.0 / (1.0 + exp(a_coeff * (kt - b_coeff)))

  # Calculate DHI
  dhi <- df * GHI

  # Calculate DNI
  cos_zenith <- cos(theta_z_deg * pi / 180)
  dni <- (GHI - dhi) / cos_zenith

  # Handle bad values
  bad_values <- (theta_z_deg > max_zenith) | (GHI < 0) | (dni < 0)
  dni <- ifelse(bad_values, 0, dni)
  # Ensure closure relationship remains valid
  dhi <- ifelse(bad_values, GHI, dhi)

  data.frame(
    time = time,
    GHI = GHI,
    DNI = dni,
    DHI = dhi,
    kt = kt,
    zenith = theta_z_deg,
    df = df
  )
}
