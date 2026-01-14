#' @title PV DC Power Pipeline: Olmo + Skoplaki + PVWatts
#'
#' @description Convenience function that computes DC power by chaining together:
#' \enumerate{
#'   \item \strong{Olmo et al. transposition}: GHI to POA irradiance
#'   \item \strong{Skoplaki cell temperature}: POA + weather to cell temperature
#'   \item \strong{PVWatts DC model}: POA + cell temp to DC power (with optional IAM)
#' }
#'
#' This function provides a complete DC power calculation pipeline. For more
#' control over individual steps, use the underlying functions directly:
#' \code{\link{olmo_transposition}}, \code{\link{skoplaki_cell_temperature}},
#' and \code{\link{pvwatts_dc}}.
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
#' @param iam_exp IAM exponent for power-law model. Default 0.05. Set to NA or
#'   FALSE to disable IAM correction.
#' @param P_dc0 DC nameplate power (W, default 230 for Trina TSM-230 PC05 module)
#' @param gamma Temperature coefficient of max power (1/K, default -0.0043)
#' @param skoplaki_variant Either "model1" or "model2" (default "model1")
#' @param T_NOCT Nominal Operating Cell Temperature (deg C, default 45)
#' @param T_a_NOCT Ambient temperature at NOCT conditions (deg C, default 20)
#' @param I_NOCT Irradiance at NOCT conditions (W/m^2, default 800)
#' @param v_NOCT Wind speed at NOCT conditions (m/s, default 1)
#' @param eta_STC Module efficiency at STC (default 0.141)
#' @param tau_alpha Product of transmittance and absorption coefficient (default 0.9)
#'
#' @return Data frame with columns: time, GHI, G_poa, T_air, wind, T_cell, P_dc,
#' zenith, sun_azimuth, incidence, skoplaki, iam (if IAM enabled)
#'
#' @seealso
#' \code{\link{olmo_transposition}} for transposition model details
#' \code{\link{skoplaki_cell_temperature}} for cell temperature model details
#' \code{\link{pvwatts_dc}} for DC power model details
#' \code{\link{pv_power_pipeline}} for complete DC + AC pipeline
#'
#' @references
#' Olmo, F. J., Vida, J., Foyo, I., Castro-Diez, Y., & Alados-Arboledas, L. (1999).
#' Prediction of global irradiance on inclined surfaces from horizontal global irradiance.
#' Energy, 24(8), 689-704. \doi{10.1016/S0360-5442(99)00025-0}
#'
#' Ayvazoğluyüksel, Ö., & Başaran Filik, Ü. (2018). Estimation methods of global
#' solar radiation, cell temperature and solar power forecasting: A review and
#' case study in Eskişehir. Renewable and Sustainable Energy Reviews, 91, 639-653.
#' \doi{10.1016/j.rser.2018.03.084}
#'
#' Skoplaki, E., Boudouvis, A. G., & Palyvos, J. A. (2008). A simple correlation
#' for the operating temperature of photovoltaic modules of arbitrary mounting.
#' Solar Energy Materials and Solar Cells, 92(11), 1393-1402.
#' \doi{10.1016/j.solmat.2008.05.016}
#'
#' Martin, N., & Ruiz, J. M. (2001). Calculation of the PV modules angular
#' losses under field conditions by means of an analytical model. Solar Energy
#' Materials and Solar Cells, 70(1), 25-38. \doi{10.1016/S0927-0248(00)00404-5}
#'
#' @export
#'
pv_dc_olmo_skoplaki_pvwatts <- function(
  time,
  lat, lon,
  GHI,
  T_air,
  wind,
  tilt,
  azimuth,
  albedo = 0.2,
  iam_exp = 0.05,
  P_dc0 = 230,
  gamma = -0.0043,
  skoplaki_variant = c("model1", "model2"),
  T_NOCT = 45,
  T_a_NOCT = 20,
  I_NOCT = 800,
  v_NOCT = 1,
  eta_STC = 0.141,
  tau_alpha = 0.9
) {
  skoplaki_variant <- match.arg(skoplaki_variant)
  stopifnot(
    length(time) == length(GHI),
    length(GHI) == length(T_air),
    length(T_air) == length(wind)
  )

  # Step 1: Olmo transposition (GHI -> G_poa)
  olmo_out <- olmo_transposition(
    time = time,
    lat = lat,
    lon = lon,
    GHI = GHI,
    tilt = tilt,
    azimuth = azimuth,
    albedo = albedo
  )

  # Step 2: Skoplaki cell temperature (G_poa, T_air, wind -> T_cell)
  T_cell <- skoplaki_cell_temperature(
    G_poa = olmo_out$G_poa,
    T_air = T_air,
    wind = wind,
    variant = skoplaki_variant,
    gamma = gamma,
    T_NOCT = T_NOCT,
    T_a_NOCT = T_a_NOCT,
    I_NOCT = I_NOCT,
    v_NOCT = v_NOCT,
    eta_STC = eta_STC,
    tau_alpha = tau_alpha
  )

  # Step 3: PVWatts DC power (G_poa, T_cell, incidence -> P_dc)
  # Only apply IAM if iam_exp is not NA/FALSE
  if (isFALSE(iam_exp) || is.na(iam_exp)) {
    incidence_param <- NULL
    iam_exp_param <- NA
    iam_values <- NA
  } else {
    incidence_param <- olmo_out$incidence
    iam_exp_param <- iam_exp
    # Calculate IAM for output
    cos_theta <- pmax(0, cos(olmo_out$incidence * pi / 180))
    iam_values <- cos_theta ^ iam_exp
  }

  P_dc <- pvwatts_dc(
    G_poa = olmo_out$G_poa,
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
    G_poa = olmo_out$G_poa,
    T_air = T_air,
    wind = wind,
    T_cell = T_cell,
    P_dc = P_dc,
    zenith = olmo_out$zenith,
    sun_azimuth = olmo_out$sun_azimuth,
    incidence = olmo_out$incidence,
    skoplaki = skoplaki_variant
  )

  # Add IAM column if IAM was enabled
  if (!isFALSE(iam_exp) && !is.na(iam_exp)) {
    result$iam <- iam_values
  }

  result
}
