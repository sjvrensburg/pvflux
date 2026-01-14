#' @title Olmo et al. Transposition Model
#'
#' @description Converts global horizontal irradiance (GHI) to plane-of-array (POA)
#' irradiance using the Olmo et al. (1999) clearness index method.
#'
#' Unlike traditional transposition models that require decomposition of GHI into
#' direct and diffuse components, this model uses the clearness index to directly
#' estimate POA irradiance. Following Ayvazoğluyüksel & Başaran Filik (2018),
#' the model is implemented as:
#' \deqn{I_{\gamma} = I \cdot \psi_o \cdot F_c}
#' where:
#' \itemize{
#'   \item \eqn{\psi_o = \exp\left(-k_t\left(\left(\frac{\pi\theta}{180}\right)^2 - \left(\frac{\pi\theta_z}{180}\right)^2\right)\right)} is the conversion function
#'   \item \eqn{F_c = 1 + \rho \sin^2(\theta/2)} is the ground reflection multiplying factor
#' }
#'
#' Note: This implementation uses angles in radians for the calculation, which is
#' mathematically equivalent to the paper's formulation using angles in degrees.
#'
#' @param time POSIXct vector of times (UTC recommended)
#' @param lat Latitude in degrees
#' @param lon Longitude in degrees
#' @param GHI Global horizontal irradiance (W/m^2)
#' @param tilt Panel tilt angle (degrees)
#' @param azimuth Panel azimuth (degrees, 0 = north)
#' @param albedo Ground albedo (default 0.2). Note: Olmo et al. used 0.25 in their
#'   original study at Granada, Spain.
#'
#' @return Data frame with columns:
#' \describe{
#'   \item{time}{Input timestamps}
#'   \item{GHI}{Input global horizontal irradiance (W/m^2)}
#'   \item{G_poa}{Plane-of-array irradiance (W/m^2)}
#'   \item{zenith}{Solar zenith angle (degrees)}
#'   \item{sun_azimuth}{Solar azimuth angle (degrees)}
#'   \item{incidence}{Angle of incidence on panel (degrees)}
#'   \item{k_t}{Clearness index}
#'   \item{I_0}{Extraterrestrial irradiance (W/m^2)}
#' }
#'
#' @references
#' Olmo, F. J., Vida, J., Foyo, I., Castro-Diez, Y., & Alados-Arboledas, L. (1999).
#' Prediction of global irradiance on inclined surfaces from horizontal global irradiance.
#' Energy, 24(8), 689-704. \doi{10.1016/S0360-5442(99)00025-0}
#'
#' @examples
#' \dontrun{
#' time <- seq(as.POSIXct("2026-01-15 06:00", tz = "UTC"),
#'             by = "hour", length.out = 12)
#' GHI <- c(50, 200, 450, 700, 850, 950, 1000, 950, 850, 700, 450, 200)
#'
#' result <- olmo_transposition(
#'   time = time,
#'   lat = -30.6279,
#'   lon = 24.0054,
#'   GHI = GHI,
#'   tilt = 20,
#'   azimuth = 0
#' )
#' }
#'
#' @export
#'
olmo_transposition <- function(
  time,
  lat,
  lon,
  GHI,
  tilt,
  azimuth,
  albedo = 0.2
) {
  stopifnot(length(time) == length(GHI))


  # Convert time to Julian Day
  jd <- JD(time)

  # Extract timezone from POSIXct object (in hours)
  tz_offset <- as.numeric(format(time[1], "%z")) / 100

  # Calculate sun vector and position
  sv <- sunvector(jd, lat, lon, tz_offset)
  sp <- sunpos(sv)

  # Extract zenith and azimuth from sunpos (column 2 = zenith, column 1 = azimuth)
  theta_z_deg <- sp[, 2]  # Solar zenith angle in degrees
  sun_az_deg <- sp[, 1]   # Solar azimuth angle in degrees
  theta_z <- theta_z_deg * pi / 180  # Convert to radians

  # Calculate panel normal vector
  beta <- tilt * pi / 180
  az <- azimuth * pi / 180
  n <- c(
    sin(beta) * sin(az),
    sin(beta) * cos(az),
    cos(beta)
  )

  # Calculate angle of incidence (theta)
  theta <- acos(pmax(-1, pmin(1, sv %*% n)))  # Incidence angle in radians
  theta_deg <- theta * 180 / pi  # Convert to degrees

  # =========================================================================
  # Olmo et al. (1999) Transposition Model - Equations 33-39
  # Reference: Olmo et al., Energy 24(8):689-704, 1999
  # =========================================================================

  # Solar constant (W/m^2)
  I_sol <- 1367

  # Day of year
  doy <- as.numeric(format(time, "%j"))

  # Eccentricity correction factor (Eq. 33 component)
  eccentricity <- 1 + 0.033 * cos(2 * pi * doy / 365)

  # Cosine of zenith angle

  cosz <- cos(theta_z)

  # Equation 33: Extraterrestrial radiation on horizontal surface
  # Paper formula: I_0 = I_sol * (1 + 0.033*cos(360d/365)) * (cos(phi)*cos(delta)*cos(W) + sin(phi)*sin(delta))
  # The trigonometric term is the standard definition for cos(theta_z). We use insol's
  # Meeus-based zenith angle calculation for higher accuracy.
  I_0 <- I_sol * eccentricity * pmax(0, cosz)

  # Equation 36: Clearness index
  k_t <- ifelse(I_0 > 0, pmin(GHI / I_0, 1), 0)

  # Equation 37: Multiplying factor (ground reflection)
  F_c <- 1 + albedo * sin(theta / 2)^2

  # Equation 38: Conversion function psi_o
  # Paper: psi_o = exp(-k_t * ((pi*theta/180)^2 - (pi*theta_z/180)^2)) where theta, theta_z are in degrees
  # Implementation: theta, theta_z are already in radians, so the formula is equivalent
  psi_o <- exp(-k_t * (theta^2 - theta_z^2))

  # Equation 39: Global radiation on inclined surface
  G_poa <- pmax(0, GHI * psi_o * F_c)

  # Handle nighttime (when sun is below horizon)
  G_poa[cosz <= 0] <- 0

  data.frame(
    time = time,
    GHI = GHI,
    G_poa = G_poa,
    zenith = theta_z_deg,
    sun_azimuth = sun_az_deg,
    incidence = theta_deg,
    k_t = k_t,
    I_0 = I_0
  )
}
