#' @title Skoplaki Cell Temperature Model
#'
#' @description Estimates PV cell temperature using the Skoplaki et al. model
#' (Equation 41 from Ayvazoğluyüksel & Başaran Filik, 2018).
#'
#' This model accounts for plane-of-array irradiance, ambient temperature,
#' wind speed effects on convective cooling, and module optical/thermal properties.
#'
#' Two wind convection coefficient variants are available:
#' model1 uses h_w = 8.91 + 2.00 v_f (Equation 42), while
#' model2 uses h_w = 5.7 + 3.8 v_w where v_w = 0.68 v_f - 0.5 (Equations 43-44).
#'
#' @param G_poa Plane-of-array irradiance (W/m^2)
#' @param T_air Ambient air temperature (deg C)
#' @param wind Wind speed measured 10m above ground (m/s)
#' @param variant Either "model1" or "model2" (default "model1")
#' @param gamma Temperature coefficient of max power (1/K, default -0.0043)
#' @param T_NOCT Nominal Operating Cell Temperature (deg C, default 45)
#' @param T_a_NOCT Ambient temperature at NOCT conditions (deg C, default 20)
#' @param I_NOCT Irradiance at NOCT conditions (W/m^2, default 800)
#' @param v_NOCT Wind speed at NOCT conditions (m/s, default 1)
#' @param eta_STC Module efficiency at STC (default 0.141)
#' @param tau_alpha Product of transmittance and absorption coefficient (default 0.9)
#'
#' @return Numeric vector of cell temperatures (deg C)
#'
#' @references
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
#' @examples
#' \dontrun{
#' G_poa <- c(200, 500, 800, 1000)
#' T_air <- c(20, 25, 30, 32)
#' wind <- c(2, 3, 4, 5)
#'
#' T_cell <- skoplaki_cell_temperature(G_poa, T_air, wind)
#' }
#'
#' @export
#'
skoplaki_cell_temperature <- function(
  G_poa,
  T_air,
  wind,
  variant = c("model1", "model2"),
  gamma = -0.0043,
  T_NOCT = 45,
  T_a_NOCT = 20,
  I_NOCT = 800,
  v_NOCT = 1,
  eta_STC = 0.141,
  tau_alpha = 0.9
) {
  variant <- match.arg(variant)
  stopifnot(
    length(G_poa) == length(T_air),
    length(T_air) == length(wind)
  )

  # Calculate wind convection coefficient at NOCT and operating conditions
  if (variant == "model1") {
    # Equation 42: h_w = 8.91 + 2.00 * v_f
    h_w_NOCT <- 8.91 + 2.00 * v_NOCT
    h_w <- 8.91 + 2.00 * wind
  } else {  # model2
    # Equations 43-44: v_w = 0.68 * v_f - 0.5; h_w = 5.7 + 3.8 * v_w
    v_w_NOCT <- pmax(0, 0.68 * v_NOCT - 0.5)
    h_w_NOCT <- 5.7 + 3.8 * v_w_NOCT
    v_w <- pmax(0, 0.68 * wind - 0.5)
    h_w <- 5.7 + 3.8 * v_w
  }

  # Skoplaki cell temperature model (Equation 41)
  T_STC <- 25  # Standard test condition temperature

  # Intermediate calculations
  irrad_ratio <- G_poa / I_NOCT
  temp_diff <- T_NOCT - T_a_NOCT
  wind_ratio <- h_w_NOCT / h_w

  # Numerator term: (1 - (eta_STC/tau_alpha) * (1 - gamma * T_STC))
  numerator_factor <- 1 - (eta_STC / tau_alpha) * (1 - gamma * T_STC)

  numerator <- T_air + irrad_ratio * temp_diff * wind_ratio * numerator_factor

  # Denominator: 1 - (gamma * eta_STC / tau_alpha) * ratio_terms
  denominator <- 1 - (gamma * eta_STC / tau_alpha) * irrad_ratio * temp_diff * wind_ratio

  T_cell <- numerator / denominator

  return(T_cell)
}
