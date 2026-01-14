#' @title PVWatts DC Power Model
#'
#' @description Calculates DC power output using the PVWatts model with
#' temperature correction.
#'
#' The model computes DC power as:
#' \deqn{P_{dc} = P_{dc0} \times \frac{G_{poa}}{1000} \times [1 + \gamma (T_{cell} - 25)]}
#'
#' where \eqn{P_{dc0}} is the nameplate DC power at STC (1000 W/m^2, 25 deg C),
#' \eqn{G_{poa}} is the plane-of-array irradiance, \eqn{\gamma} is the temperature
#' coefficient, and \eqn{T_{cell}} is the cell temperature.
#'
#' @param G_poa Plane-of-array irradiance (W/m^2)
#' @param T_cell Cell temperature (deg C)
#' @param P_dc0 DC nameplate power at STC (W, default 230 for Trina TSM-230 PC05)
#' @param gamma Temperature coefficient of max power (1/K, default -0.0043)
#'
#' @return Numeric vector of DC power output (W). Negative values are set to zero.
#'
#' @references
#' Dobos, A. P. (2014). PVWatts Version 5 Manual. NREL/TP-6A20-62641.
#' National Renewable Energy Laboratory.
#'
#' @examples
#' \dontrun{
#' G_poa <- c(200, 500, 800, 1000)
#' T_cell <- c(30, 40, 50, 55)
#'
#' # Single module
#' P_dc <- pvwatts_dc(G_poa, T_cell, P_dc0 = 230)
#'
#' # Full plant (44,880 modules)
#' P_dc_plant <- pvwatts_dc(G_poa, T_cell, P_dc0 = 44880 * 230)
#' }
#'
#' @export
#'
pvwatts_dc <- function(
  G_poa,
  T_cell,
  P_dc0 = 230,
  gamma = -0.0043
) {
  stopifnot(length(G_poa) == length(T_cell))

  # PVWatts DC power equation
  P_dc <- P_dc0 * (G_poa / 1000) * (1 + gamma * (T_cell - 25))

  # Ensure non-negative output

  P_dc[P_dc < 0] <- 0

  return(P_dc)
}
