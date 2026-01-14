#' @title PV DC Power Pipeline: Olmo + Skoplaki + PVWatts
#'
#' @description Computes DC power using Olmo transposition for POA irradiance, Skoplaki cell temperature, and PVWatts model.
#'
#' @param time POSIXct vector of times (UTC recommended)
#'
#' @param lat Latitude in degrees
#'
#' @param lon Longitude in degrees
#'
#' @param GHI Global horizontal irradiance (W/m^2)
#'
#' @param T_air Ambient air temperature (°C)
#'
#' @param wind Wind speed (m/s)
#'
#' @param tilt Panel tilt angle (degrees)
#'
#' @param azimuth Panel azimuth (degrees, 0 = north)
#'
#' @param albedo Ground albedo (default 0.2)
#'
#' @param P_dc0 DC nameplate power (W, default 230 for Trina TSM-230 PC05 module)
#'
#' @param gamma Temperature coefficient of max power (1/K, default -0.0043 for TSM-230)
#'
#' @param skoplaki_variant Either "model1" or "model2" (default "model1"). Model 1 uses h_w = 8.91 + 2.00*v_f, Model 2 uses h_w = 5.7 + 3.8*v_w
#'
#' @param T_NOCT Nominal Operating Cell Temperature in °C (default 45 for TSM-230)
#'
#' @param T_a_NOCT Ambient temperature at NOCT conditions in °C (default 20)
#'
#' @param I_NOCT Irradiance at NOCT conditions in W/m² (default 800)
#'
#' @param v_NOCT Wind speed at NOCT conditions in m/s (default 1)
#'
#' @param eta_STC Module efficiency at STC (default 0.141 for TSM-230)
#'
#' @param tau_alpha Product of transmittance and absorption coefficient (default 0.9)
#'
#' @return Data frame with G_poa, T_cell, P_dc, etc.
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
  # Convert time to Julian Day
  jd <- JD(time)
  # Extract timezone from POSIXct object (in hours)
  # If time is UTC, tz_offset will be 0
  tz_offset <- as.numeric(format(time[1], "%z")) / 100
  # Calculate sun vector and position
  sv <- sunvector(jd, lat, lon, tz_offset)
  sp <- sunpos(sv)
  # Extract zenith from matrix (column 2)
  theta_z <- sp[, 2] * pi/180
  beta <- tilt * pi/180
  az <- azimuth * pi/180
  n <- c(
    sin(beta) * sin(az),
    sin(beta) * cos(az),
    cos(beta)
  )
  theta <- acos(pmax(-1, pmin(1, sv %*% n)))
  cosz <- cos(theta_z)
  cost <- cos(theta)
  R_T <- ifelse(
    cosz > 0,
    cost / cosz + albedo * (1 - cos(beta)) / 2,
    0
  )
  G_poa <- pmax(0, GHI * R_T)

  # Calculate wind convection coefficient at NOCT
  if (skoplaki_variant == "model1") {
    h_w_NOCT <- 8.91 + 2.00 * v_NOCT
    h_w <- 8.91 + 2.00 * wind
  } else {  # model2
    v_w_NOCT <- 0.68 * v_NOCT - 0.5
    h_w_NOCT <- 5.7 + 3.8 * v_w_NOCT
    v_w <- 0.68 * wind - 0.5
    h_w <- 5.7 + 3.8 * v_w
  }

  # Skoplaki cell temperature model (Equation 41 from paper)
  # Reference: Skoplaki et al. (2008), DOI: 10.1016/j.solmat.2008.05.016
  T_STC <- 25  # Standard test condition temperature

  # Numerator
  numerator <- T_air + (G_poa / I_NOCT) * (T_NOCT - T_a_NOCT) *
    (h_w_NOCT / h_w) * (1 - eta_STC * (1 - gamma * T_STC)) / tau_alpha

  # Denominator
  denominator <- 1 - (gamma * eta_STC / tau_alpha) *
    (G_poa / I_NOCT) * (T_NOCT - T_a_NOCT) * (h_w_NOCT / h_w)

  T_cell <- numerator / denominator

  P_dc <- P_dc0 * (G_poa / 1000) * (1 + gamma * (T_cell - 25))
  P_dc[P_dc < 0] <- 0
  data.frame(
    time = time,
    GHI = GHI,
    G_poa = G_poa,
    T_air = T_air,
    wind = wind,
    T_cell = T_cell,
    P_dc = P_dc,
    zenith = sp[, 2],
    azimuth = sp[, 1],
    incidence = theta * 180/pi,
    skoplaki = skoplaki_variant
  )
}
