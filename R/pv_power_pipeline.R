#' @title Complete PV Power Pipeline: DC and AC
#'
#' @description Convenience function that calculates both DC and AC power by
#' chaining together all four models:
#' \enumerate{
#'   \item \strong{Olmo et al. transposition}: GHI to POA irradiance
#'   \item \strong{Skoplaki cell temperature}: POA + weather to cell temperature
#'   \item \strong{PVWatts DC model}: POA + cell temp to DC power
#'   \item \strong{Simple AC clipping}: DC power to AC power with inverter losses
#' }
#'
#' This is the highest-level convenience function in the package. For more control
#' over individual steps, use the underlying functions:
#' \itemize{
#'   \item \code{\link{olmo_transposition}} - Transposition model only
#'   \item \code{\link{skoplaki_cell_temperature}} - Cell temperature only
#'   \item \code{\link{pvwatts_dc}} - DC power calculation only
#'   \item \code{\link{pv_ac_simple_clipping}} - AC conversion only
#'   \item \code{\link{pv_dc_olmo_skoplaki_pvwatts}} - DC pipeline (steps 1-3)
#' }
#'
#' @param time POSIXct vector of times (UTC recommended)
#' @param lat Latitude in degrees
#' @param lon Longitude in degrees
#' @param GHI Global horizontal irradiance (W/m^2)
#' @param T_air Ambient air temperature (deg C)
#' @param wind Wind speed (m/s)
#' @param tilt Panel tilt angle (degrees)
#' @param azimuth Panel azimuth (degrees, 0 = north)
#' @param albedo Ground albedo (default 0.2)
#' @param P_dc0 DC nameplate power (W, default 230 for Trina TSM-230 PC05 module)
#' @param gamma Temperature coefficient of max power (1/K, default -0.0043)
#' @param skoplaki_variant Either "model1" or "model2" (default "model1")
#' @param T_NOCT Nominal Operating Cell Temperature (deg C, default 45)
#' @param T_a_NOCT Ambient temperature at NOCT conditions (deg C, default 20)
#' @param I_NOCT Irradiance at NOCT conditions (W/m^2, default 800)
#' @param v_NOCT Wind speed at NOCT conditions (m/s, default 1)
#' @param eta_STC Module efficiency at STC (default 0.141)
#' @param tau_alpha Product of transmittance and absorption coefficient (default 0.9)
#' @param n_inverters Number of inverters (default 20)
#' @param inverter_kw kW rating per inverter (default 500)
#' @param eta_inv Inverter efficiency (default 0.97)
#'
#' @return Data frame with columns: time, GHI, G_poa, T_air, wind, T_cell, P_dc,
#' P_ac, clipped, P_ac_rated, zenith, sun_azimuth, incidence, skoplaki
#'
#' @seealso
#' \code{\link{olmo_transposition}} for transposition model details
#' \code{\link{skoplaki_cell_temperature}} for cell temperature model details
#' \code{\link{pvwatts_dc}} for DC power model details
#' \code{\link{pv_ac_simple_clipping}} for AC conversion details
#' \code{\link{pv_dc_olmo_skoplaki_pvwatts}} for DC-only pipeline
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example for De Aar plant (10.32 MW)
#' time <- seq(as.POSIXct("2026-01-15 06:00", tz = "UTC"),
#'             by = "hour", length.out = 12)
#' GHI <- c(50, 200, 450, 700, 850, 950, 1000, 950, 850, 700, 450, 200)
#' T_air <- c(20, 22, 26, 29, 32, 34, 35, 34, 32, 29, 26, 23)
#' wind <- c(2, 2.5, 3, 4, 4.5, 5, 5, 5, 4.5, 4, 3, 2.5)
#'
#' result <- pv_power_pipeline(
#'   time = time,
#'   lat = -30.6279,
#'   lon = 24.0054,
#'   GHI = GHI,
#'   T_air = T_air,
#'   wind = wind,
#'   tilt = 20,
#'   azimuth = 0,
#'   P_dc0 = 44880 * 230,  # 44,880 modules x 230W
#'   n_inverters = 20,
#'   inverter_kw = 500
#' )
#'
#' head(result)
#' }
#'
pv_power_pipeline <- function(
  time,
  lat, lon,
  GHI,
  T_air,
  wind,
  tilt,
  azimuth,
  albedo = 0.2,
  P_dc0 = 230,
  gamma = -0.0043,
  skoplaki_variant = c("model1", "model2"),
  T_NOCT = 45,
  T_a_NOCT = 20,
  I_NOCT = 800,
  v_NOCT = 1,
  eta_STC = 0.141,
  tau_alpha = 0.9,
  n_inverters = 20,
  inverter_kw = 500,
  eta_inv = 0.97
) {
  # Calculate DC power using the DC pipeline
  dc_out <- pv_dc_olmo_skoplaki_pvwatts(
    time = time,
    lat = lat,
    lon = lon,
    GHI = GHI,
    T_air = T_air,
    wind = wind,
    tilt = tilt,
    azimuth = azimuth,
    albedo = albedo,
    P_dc0 = P_dc0,
    gamma = gamma,
    skoplaki_variant = skoplaki_variant,
    T_NOCT = T_NOCT,
    T_a_NOCT = T_a_NOCT,
    I_NOCT = I_NOCT,
    v_NOCT = v_NOCT,
    eta_STC = eta_STC,
    tau_alpha = tau_alpha
  )

  # Calculate AC power using simple clipping model
  ac_out <- pv_ac_simple_clipping(
    P_dc = dc_out$P_dc,
    n_inverters = n_inverters,
    inverter_kw = inverter_kw,
    eta_inv = eta_inv
  )

  # Combine results
  dc_out$P_ac <- ac_out$P_ac
  dc_out$clipped <- ac_out$clipped
  dc_out$P_ac_rated <- ac_out$P_ac_rated

  return(dc_out)
}
