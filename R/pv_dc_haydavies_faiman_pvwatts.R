#' @title PV DC Power Pipeline: Erbs + Hay-Davies + Faiman + PVWatts
#'
#' @description Convenience function that computes DC power by chaining together
#' the alternative model set:
#' \enumerate{
#'   \item \strong{Erbs decomposition}: GHI to DNI and DHI
#'   \item \strong{Hay-Davies transposition}: GHI, DNI, DHI to POA irradiance
#'   \item \strong{Faiman cell temperature}: POA + weather to cell temperature
#'   \item \strong{PVWatts DC model}: POA + cell temp to DC power (with optional IAM)
#' }
#'
#' This function provides an alternative DC power calculation pipeline using
#' well-established models from the pvlib library. For more control over
#' individual steps, use the underlying functions directly:
#' \code{\link{erbs_decomposition}}, \code{\link{haydavies_transposition}},
#' \code{\link{faiman_cell_temperature}}, and \code{\link{pvwatts_dc}}.
#'
#' @param time Timestamps as POSIXct, POSIXlt, character, or numeric. If a timezone
#'   is specified, times are internally converted to UTC for solar position
#'   calculations and returned in the original timezone. If no timezone is
#'   specified, UTC is assumed. See \code{\link{time_utils}} for details.
#' @param lat Latitude in degrees
#' @param lon Longitude in degrees
#' @param GHI Global horizontal irradiance (W/m^2)
#' @param T_air Ambient air temperature (deg C)
#' @param wind Wind speed (m/s)
#' @param tilt Panel tilt angle (degrees)
#' @param azimuth Panel azimuth (degrees, 0 = north)
#' @param albedo Ground albedo (default 0.2)
#' @param iam_exp IAM exponent for power-law model. Default 0.05. Set to NA or
#'   FALSE to disable IAM correction.
#' @param u0 Combined heat loss factor coefficient for Faiman model.
#'   Default 25.0 W/(m²·°C).
#' @param u1 Combined heat loss factor influenced by wind for Faiman model.
#'   Default 6.84 W/(m²·°C·m/s).
#' @param min_cos_zenith Minimum value of cos(zenith) for Hay-Davies Rb calculation.
#'   Default: 0.01745
#' @param P_dc0 DC nameplate power (W, default 230 for Trina TSM-230 PC05 module)
#' @param gamma Temperature coefficient of max power (1/K, default -0.0043)
#'
#' @return Data frame with columns: time, GHI, DNI, DHI, G_poa, T_air, wind,
#' T_cell, P_dc, zenith, azimuth, incidence, ai, rb, iam (if IAM enabled)
#'
#' @seealso
#' \code{\link{erbs_decomposition}} for Erbs decomposition model details
#' \code{\link{haydavies_transposition}} for Hay-Davies transposition details
#' \code{\link{faiman_cell_temperature}} for Faiman cell temperature details
#' \code{\link{pvwatts_dc}} for DC power model details
#' \code{\link{pv_dc_olmo_skoplaki_pvwatts}} for Olmo/Skoplaki pipeline
#' \code{\link{pv_power_pipeline}} for complete DC + AC pipeline
#'
#' @references
#' Erbs, D. G., Klein, S. A., & Duffie, J. A. (1982). Estimation of the
#' diffuse radiation fraction for hourly, daily and monthly-average global
#' radiation. Solar Energy, 28(4), 293-302.
#'
#' Hay, J. E., & Davies, J. A. (1980). Calculations of the solar radiation
#' incident on an inclined surface. In J. E. Hay & T. K. Won (Eds.),
#' Proc. of First Canadian Solar Radiation Data Workshop (pp. 59).
#' Ministry of Supply and Services, Canada.
#'
#' Faiman, D. (2008). Assessing the outdoor operating temperature of
#' photovoltaic modules. Progress in Photovoltaics, 16(4), 307-315.
#'
#' @examples
#' \dontrun{
#' time <- seq(as.POSIXct("2026-01-15 08:00", tz = "UTC"),
#'             by = "hour", length.out = 6)
#' GHI <- c(450, 700, 850, 950, 850, 700)
#' T_air <- c(26, 29, 32, 34, 32, 29)
#' wind <- c(3, 4, 4.5, 5, 4.5, 4)
#'
#' result <- pv_dc_haydavies_faiman_pvwatts(
#'   time = time,
#'   lat = -30.6279,
#'   lon = 24.0054,
#'   GHI = GHI,
#'   T_air = T_air,
#'   wind = wind,
#'   tilt = 20,
#'   azimuth = 0
#' )
#'
#' head(result)
#' }
#'
#' @export
#'
pv_dc_haydavies_faiman_pvwatts <- function(
  time,
  lat, lon,
  GHI,
  T_air,
  wind,
  tilt,
  azimuth,
  albedo = 0.2,
  iam_exp = 0.05,
  u0 = 25.0,
  u1 = 6.84,
  min_cos_zenith = 0.01745,
  P_dc0 = 230,
  gamma = -0.0043
) {
  stopifnot(
    length(time) == length(GHI),
    length(GHI) == length(T_air),
    length(T_air) == length(wind)
  )

  # Step 1: Erbs decomposition (GHI -> DNI, DHI)
  erbs_out <- erbs_decomposition(
    time = time,
    lat = lat,
    lon = lon,
    GHI = GHI
  )

  # Step 2: Hay-Davies transposition (GHI, DNI, DHI -> POA)
  haydavies_out <- haydavies_transposition(
    time = time,
    lat = lat,
    lon = lon,
    GHI = GHI,
    DNI = erbs_out$DNI,
    DHI = erbs_out$DHI,
    tilt = tilt,
    azimuth = azimuth,
    albedo = albedo,
    min_cos_zenith = min_cos_zenith
  )

  # Step 3: Faiman cell temperature (POA, T_air, wind -> T_cell)
  T_cell <- faiman_cell_temperature(
    poa_global = haydavies_out$poa_global,
    temp_air = T_air,
    wind_speed = wind,
    u0 = u0,
    u1 = u1
  )

  # Step 4: PVWatts DC power (POA, T_cell, incidence -> P_dc)
  # Only apply IAM if iam_exp is not NA/FALSE
  if (isFALSE(iam_exp) || is.na(iam_exp)) {
    incidence_param <- NULL
    iam_exp_param <- NA
    iam_values <- NA
  } else {
    incidence_param <- haydavies_out$incidence
    iam_exp_param <- iam_exp
    # Calculate IAM for output
    cos_theta <- pmax(0, cos(haydavies_out$incidence * pi / 180))
    iam_values <- cos_theta ^ iam_exp
  }

  P_dc <- pvwatts_dc(
    G_poa = haydavies_out$poa_global,
    T_cell = T_cell,
    incidence = incidence_param,
    iam_exp = iam_exp_param,
    P_dc0 = P_dc0,
    gamma = gamma
  )

  # Combine results
  result <- data.frame(
    time = time,
    GHI = GHI,
    DNI = erbs_out$DNI,
    DHI = erbs_out$DHI,
    G_poa = haydavies_out$poa_global,
    T_air = T_air,
    wind = wind,
    T_cell = T_cell,
    P_dc = P_dc,
    zenith = haydavies_out$zenith,
    azimuth = haydavies_out$azimuth,
    incidence = haydavies_out$incidence,
    ai = haydavies_out$ai,
    rb = haydavies_out$rb
  )

  # Add IAM column if IAM was enabled
  if (!isFALSE(iam_exp) && !is.na(iam_exp)) {
    result$iam <- iam_values
  }

  result
}
