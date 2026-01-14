#' Faiman Cell Temperature Model
#'
#' @description
#' Calculate PV cell or module temperature using the Faiman (2008) model.
#'
#' The Faiman model uses an empirical heat loss factor model and is adopted in
#' the IEC 61853 standards. The model calculates the cell temperature as:
#'
#' \deqn{T_{cell} = T_{amb} + \frac{G_{POA}}{u_0 + u_1 \times v}}
#'
#' where \eqn{G_{POA}} is the plane-of-array irradiance, \eqn{v} is the wind
#' speed, and \eqn{u_0} and \eqn{u_1} are the combined heat loss factor
#' coefficients.
#'
#' @param poa_global Plane-of-array irradiance (W/m^2)
#' @param temp_air Ambient dry bulb temperature (°C)
#' @param wind_speed Wind speed in m/s measured at the same height for which
#'   the wind loss factor was determined. Default: 1.0 (m/s at module height
#'   used to determine NOCT)
#' @param u0 Combined heat loss factor coefficient. Default: 25.0 W/(m²·°C),
#'   determined by Faiman for 7 silicon modules in the Negev desert on an open
#'   rack at 30.9° tilt.
#' @param u1 Combined heat loss factor influenced by wind. Default: 6.84
#'   W/(m²·°C·m/s), determined by Faiman for the same setup.
#'
#' @return Numeric vector of cell/module temperatures in degrees Celsius
#'
#' @references
#' Faiman, D. (2008). Assessing the outdoor operating temperature of
#' photovoltaic modules. Progress in Photovoltaics, 16(4), 307-315.
#' \doi{10.1002/pip.813}
#'
#' IEC 61853-2. Photovoltaic (PV) module performance testing and energy
#' rating - Part 2: Spectral responsivity, incidence angle and module
#' operating temperature measurements. IEC, Geneva, 2018.
#'
#' IEC 61853-3. Photovoltaic (PV) module performance testing and energy
#' rating - Part 3: Energy rating of PV modules. IEC, Geneva, 2018.
#'
#' @examples
#' \dontrun{
#' poa_global <- c(0, 100, 500, 800, 1000)
#' temp_air <- c(20, 22, 25, 28, 30)
#' wind_speed <- c(1.0, 1.5, 2.0, 1.0, 0.5)
#'
#' t_cell <- faiman_cell_temperature(
#'   poa_global = poa_global,
#'   temp_air = temp_air,
#'   wind_speed = wind_speed
#' )
#'
#' print(t_cell)
#' }
#'
#' @export
#'
faiman_cell_temperature <- function(
  poa_global,
  temp_air,
  wind_speed = 1.0,
  u0 = 25.0,
  u1 = 6.84
) {
  # Calculate total heat loss factor
  total_loss_factor <- u0 + u1 * wind_speed

  # Heat input is the POA irradiance
  heat_input <- poa_global

  # Temperature rise due to solar irradiance
  temp_difference <- heat_input / total_loss_factor

  # Cell temperature is ambient plus temperature rise
  temp_cell <- temp_air + temp_difference

  return(temp_cell)
}
