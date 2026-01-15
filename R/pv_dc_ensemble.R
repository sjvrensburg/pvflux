#' @title PV DC Power Ensemble - All Model Combinations
#'
#' @description Computes DC power using all combinations of transposition,
#' decomposition, and cell temperature models, returning results in a format
#' suitable for ensemble analysis.
#'
#' This function runs a complete ensemble of model combinations covering:
#' \enumerate{
#'   \item Hay-Davies + ERBS decomposition + Skoplaki (model1/model2)
#'   \item Hay-Davies + ERBS decomposition + Faiman
#'   \item Hay-Davies + Boland decomposition + Skoplaki (model1/model2)
#'   \item Hay-Davies + Boland decomposition + Faiman
#'   \item Reindl + ERBS decomposition + Skoplaki (model1/model2)
#'   \item Reindl + ERBS decomposition + Faiman
#'   \item Reindl + Boland decomposition + Skoplaki (model1/model2)
#'   \item Reindl + Boland decomposition + Faiman
#'   \item Perez + ERBS decomposition + Skoplaki (model1/model2)
#'   \item Perez + ERBS decomposition + Faiman
#'   \item Perez + Boland decomposition + Skoplaki (model1/model2)
#'   \item Perez + Boland decomposition + Faiman
#'   \item Olmo + Skoplaki (model1/model2)
#'   \item Olmo + Faiman
#' }
#'
#' The skoplaki_variant parameter only affects results when cell_temp = "skoplaki".
#' For faiman cell temperature, only one variant is run since results are identical.
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
#' @param P_dc0 DC nameplate power (W, default 230 for Trina TSM-230 PC05 module)
#' @param gamma Temperature coefficient of max power (1/K, default -0.0043)
#' @param min_cos_zenith Minimum value of cos(zenith) for Hay-Davies Rb calculation.
#'   Default: 0.01745.
#' @param T_NOCT Nominal Operating Cell Temperature (deg C, default 45).
#' @param T_a_NOCT Ambient temperature at NOCT conditions (deg C, default 20).
#' @param I_NOCT Irradiance at NOCT conditions (W/m^2, default 800).
#' @param v_NOCT Wind speed at NOCT conditions (m/s, default 1).
#' @param eta_STC Module efficiency at STC (default 0.141).
#' @param tau_alpha Product of transmittance and absorption coefficient (default 0.9).
#' @param u0 Combined heat loss factor coefficient for Faiman model.
#'   Default 25.0 W/(m²·°C).
#' @param u1 Combined heat loss factor influenced by wind for Faiman model.
#'   Default 6.84 W/(m²·°C·m/s).
#'
#' @return Data frame with columns:
#' \itemize{
#'   \item time, GHI, T_air, wind - Input data
#'   \item model - Model combination identifier (e.g., "haydavies_erbs_skoplaki_model1")
#'   \item transposition - Transposition model used
#'   \item decomposition - Decomposition model used (erbs/boland for applicable models, NA for olmo)
#'   \item cell_temp - Cell temperature model used
#'   \item skoplaki_variant - Skoplaki model variant used (model1/model2)
#'   \item G_poa - Plane-of-array irradiance (W/m^2)
#'   \item T_cell - Cell temperature (deg C)
#'   \item P_dc - DC power output (W)
#'   \item iam - Incidence angle modifier (if IAM enabled)
#' }
#'
#' @examples
#' \dontrun{
#' time <- seq(as.POSIXct("2026-01-15 08:00", tz = "UTC"),
#'             by = "hour", length.out = 6)
#' GHI <- c(450, 700, 850, 950, 850, 700)
#' T_air <- c(26, 29, 32, 34, 32, 29)
#' wind <- c(3, 4, 4.5, 5, 4.5, 4)
#'
#' ensemble_result <- pv_dc_ensemble(
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
#' # View unique model combinations
#' unique(ensemble_result$model)
#'
#' # Calculate ensemble statistics
#' ensemble_stats <- aggregate(P_dc ~ time, data = ensemble_result,
#'                            FUN = function(x) c(mean = mean(x), sd = sd(x)))
#' }
#'
#' @seealso
#' \code{\link{pv_dc_pipeline}} for single model combination
#' \code{\link{pv_power_ensemble}} for full ensemble including AC power
#'
#' @export
#'
pv_dc_ensemble <- function(
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
  min_cos_zenith = 0.01745,
  T_NOCT = 45,
  T_a_NOCT = 20,
  I_NOCT = 800,
  v_NOCT = 1,
  eta_STC = 0.141,
  tau_alpha = 0.9,
  u0 = 25.0,
  u1 = 6.84
) {
  # Define transposition models that need decomposition
  decomp_transpositions <- c("haydavies", "reindl", "perez")

  # Run ensemble: only run skoplaki_variant when it affects results
  results <- list()

  # For transposition models needing decomposition
  for (transposition in decomp_transpositions) {
    for (decomp_model in c("erbs", "boland")) {
      # Faiman cell temp (skoplaki_variant doesn't affect result)
      result <- pv_dc_pipeline(
        time = time,
        lat = lat,
        lon = lon,
        GHI = GHI,
        T_air = T_air,
        wind = wind,
        tilt = tilt,
        azimuth = azimuth,
        albedo = albedo,
        transposition_model = transposition,
        decomposition_model = decomp_model,
        cell_temp_model = "faiman",
        iam_exp = iam_exp,
        P_dc0 = P_dc0,
        gamma = gamma,
        min_cos_zenith = min_cos_zenith,
        skoplaki_variant = "model1",  # Doesn't affect faiman, but required
        T_NOCT = T_NOCT,
        T_a_NOCT = T_a_NOCT,
        I_NOCT = I_NOCT,
        v_NOCT = v_NOCT,
        eta_STC = eta_STC,
        tau_alpha = tau_alpha,
        u0 = u0,
        u1 = u1
      )
      result$model <- paste(transposition, decomp_model, "faiman", sep = "_")
      result$skoplaki_variant <- NA
      cols_to_keep <- c("time", "GHI", "T_air", "wind", "model", "transposition",
                        "decomposition", "cell_temp", "skoplaki_variant",
                        "G_poa", "T_cell", "P_dc")
      if (!isFALSE(iam_exp) && !is.na(iam_exp)) {
        cols_to_keep <- c(cols_to_keep, "iam")
      }
      results[[length(results) + 1]] <- result[, cols_to_keep]

      # Skoplaki cell temp (both variants)
      for (skoplaki_variant in c("model1", "model2")) {
        result <- pv_dc_pipeline(
          time = time,
          lat = lat,
          lon = lon,
          GHI = GHI,
          T_air = T_air,
          wind = wind,
          tilt = tilt,
          azimuth = azimuth,
          albedo = albedo,
          transposition_model = transposition,
          decomposition_model = decomp_model,
          cell_temp_model = "skoplaki",
          iam_exp = iam_exp,
          P_dc0 = P_dc0,
          gamma = gamma,
          min_cos_zenith = min_cos_zenith,
          skoplaki_variant = skoplaki_variant,
          T_NOCT = T_NOCT,
          T_a_NOCT = T_a_NOCT,
          I_NOCT = I_NOCT,
          v_NOCT = v_NOCT,
          eta_STC = eta_STC,
          tau_alpha = tau_alpha,
          u0 = u0,
          u1 = u1
        )

        # Rename skoplaki column to skoplaki_variant
        result$skoplaki_variant <- result$skoplaki
        result$skoplaki <- NULL

        result$model <- paste(transposition, decomp_model, "skoplaki", skoplaki_variant, sep = "_")

        cols_to_keep <- c("time", "GHI", "T_air", "wind", "model", "transposition",
                          "decomposition", "cell_temp", "skoplaki_variant",
                          "G_poa", "T_cell", "P_dc")
        if (!isFALSE(iam_exp) && !is.na(iam_exp)) {
          cols_to_keep <- c(cols_to_keep, "iam")
        }
        results[[length(results) + 1]] <- result[, cols_to_keep]
      }
    }
  }

  # Olmo (no decomposition needed)
  # Faiman cell temp
  result <- pv_dc_pipeline(
    time = time,
    lat = lat,
    lon = lon,
    GHI = GHI,
    T_air = T_air,
    wind = wind,
    tilt = tilt,
    azimuth = azimuth,
    albedo = albedo,
    transposition_model = "olmo",
    cell_temp_model = "faiman",
    iam_exp = iam_exp,
    P_dc0 = P_dc0,
    gamma = gamma,
    min_cos_zenith = min_cos_zenith,
    skoplaki_variant = "model1",
    T_NOCT = T_NOCT,
    T_a_NOCT = T_a_NOCT,
    I_NOCT = I_NOCT,
    v_NOCT = v_NOCT,
    eta_STC = eta_STC,
    tau_alpha = tau_alpha,
    u0 = u0,
    u1 = u1
  )
  result$model <- paste("olmo", "faiman", sep = "_")
  result$skoplaki_variant <- NA
  cols_to_keep <- c("time", "GHI", "T_air", "wind", "model", "transposition",
                    "decomposition", "cell_temp", "skoplaki_variant",
                    "G_poa", "T_cell", "P_dc")
  if (!isFALSE(iam_exp) && !is.na(iam_exp)) {
    cols_to_keep <- c(cols_to_keep, "iam")
  }
  results[[length(results) + 1]] <- result[, cols_to_keep]

  # Skoplaki cell temp (both variants)
  for (skoplaki_variant in c("model1", "model2")) {
    result <- pv_dc_pipeline(
      time = time,
      lat = lat,
      lon = lon,
      GHI = GHI,
      T_air = T_air,
      wind = wind,
      tilt = tilt,
      azimuth = azimuth,
      albedo = albedo,
      transposition_model = "olmo",
      cell_temp_model = "skoplaki",
      iam_exp = iam_exp,
      P_dc0 = P_dc0,
      gamma = gamma,
      min_cos_zenith = min_cos_zenith,
      skoplaki_variant = skoplaki_variant,
      T_NOCT = T_NOCT,
      T_a_NOCT = T_a_NOCT,
      I_NOCT = I_NOCT,
      v_NOCT = v_NOCT,
      eta_STC = eta_STC,
      tau_alpha = tau_alpha,
      u0 = u0,
      u1 = u1
    )

    # Rename skoplaki column to skoplaki_variant
    result$skoplaki_variant <- result$skoplaki
    result$skoplaki <- NULL

    result$model <- paste("olmo", "skoplaki", skoplaki_variant, sep = "_")

    cols_to_keep <- c("time", "GHI", "T_air", "wind", "model", "transposition",
                      "decomposition", "cell_temp", "skoplaki_variant",
                      "G_poa", "T_cell", "P_dc")
    if (!isFALSE(iam_exp) && !is.na(iam_exp)) {
      cols_to_keep <- c(cols_to_keep, "iam")
    }
    results[[length(results) + 1]] <- result[, cols_to_keep]
  }

  # Combine all results
  do.call(rbind, results)
}
