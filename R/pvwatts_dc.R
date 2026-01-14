#' @title PVWatts DC Power Model
#'
#' @description Calculates DC power output using the PVWatts model with
#' temperature correction and optional Incidence Angle Modifier (IAM).
#'
#' The model computes DC power as:
#' \deqn{P_{dc} = P_{dc0} \times \frac{G_{poa} \times IAM}{1000} \times [1 + \gamma (T_{cell} - 25)]}
#'
#' where \eqn{P_{dc0}} is the nameplate DC power at STC (1000 W/m^2, 25 deg C),
#' \eqn{G_{poa}} is the plane-of-array irradiance, \eqn{IAM} is the incidence
#' angle modifier, \eqn{\gamma} is the temperature coefficient, and \eqn{T_{cell}}
#' is the cell temperature.
#'
#' When enabled, the IAM uses a simple power-law model:
#' \deqn{IAM = \cos(\theta)^b}
#' where \eqn{\theta} is the incidence angle and \eqn{b} is an empirical
#' coefficient (default 0.05). This accounts for optical losses at high incidence
#' angles (>60°) and can correct 1-3% underestimation at low sun angles.
#'
#' @param G_poa Plane-of-array irradiance (W/m^2)
#' @param T_cell Cell temperature (deg C)
#' @param incidence Incidence angle on panel (degrees). Optional - if provided,
#'   IAM will be applied.
#' @param iam_exp IAM exponent for power-law model. Default 0.05. Set to NA or
#'   NULL to disable IAM correction.
#' @param P_dc0 DC nameplate power at STC (W, default 230 for Trina TSM-230 PC05)
#' @param gamma Temperature coefficient of max power (1/K, default -0.0043)
#'
#' @return Numeric vector of DC power output (W). Negative values are set to zero.
#'
#' @references
#' Dobos, A. P. (2014). PVWatts Version 5 Manual. NREL/TP-6A20-62641.
#' National Renewable Energy Laboratory.
#'
#' Martin, N., & Ruiz, J. M. (2001). Calculation of the PV modules angular
#' losses under field conditions by means of an analytical model. Solar Energy
#' Materials and Solar Cells, 70(1), 25-38. \doi{10.1016/S0927-0248(00)00404-5}
#'
#' @examples
#' \dontrun{
#' G_poa <- c(200, 500, 800, 1000)
#' T_cell <- c(30, 40, 50, 55)
#' incidence <- c(10, 25, 45, 70)
#'
#' # Without IAM
#' P_dc <- pvwatts_dc(G_poa, T_cell, P_dc0 = 230)
#'
#' # With IAM correction
#' P_dc_iam <- pvwatts_dc(G_poa, T_cell, incidence = incidence, P_dc0 = 230)
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
  incidence = NULL,
  iam_exp = 0.05,
  P_dc0 = 230,
  gamma = -0.0043
) {
  stopifnot(length(G_poa) == length(T_cell))

  # Apply IAM if incidence angle is provided and iam_exp is not NA
  if (!is.null(incidence) && !any(is.na(iam_exp))) {
    stopifnot(length(G_poa) == length(incidence))
    # Power-law IAM: cos(theta)^b
    # Clamp cos(theta) to avoid negative values for large angles
    cos_theta <- pmax(0, cos(incidence * pi / 180))
    iam <- cos_theta ^ iam_exp
    G_poa <- G_poa * iam
  }

  # PVWatts DC power equation
  P_dc <- P_dc0 * (G_poa / 1000) * (1 + gamma * (T_cell - 25))

  # Ensure non-negative output
  P_dc[P_dc < 0] <- 0

  return(P_dc)
}
