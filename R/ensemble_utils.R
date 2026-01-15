#' @title Ensemble Summary Statistics
#'
#' @description Compute summary statistics across all ensemble models for each timestep.
#'
#' @param ensemble_result Data frame returned by \code{\link{pv_dc_ensemble}} or
#'   \code{\link{pv_power_ensemble}}.
#' @param value_column Column name to compute statistics for (default: "P_dc" for DC
#'   ensemble, "P_ac" for power ensemble).
#' @param probs Numeric vector of probabilities for quantiles (default: c(0.05, 0.25, 0.5, 0.75, 0.95)).
#'
#' @return Data frame with one row per timestep and columns for each statistic.
#'
#' @examples
#' \dontrun{
#' ens <- pv_dc_ensemble(...)
#' summary <- ensemble_summary(ens, "P_dc")
#' }
#'
#' @export
#'
ensemble_summary <- function(ensemble_result, value_column = NULL, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)) {
  # Auto-detect value column if not specified
  if (is.null(value_column)) {
    if ("P_ac" %in% names(ensemble_result)) {
      value_column <- "P_ac"
    } else if ("P_dc" %in% names(ensemble_result)) {
      value_column <- "P_dc"
    } else {
      stop("Could not auto-detect value column. Please specify value_column.")
    }
  }

  if (!value_column %in% names(ensemble_result)) {
    stop(paste("Column", value_column, "not found in ensemble_result."))
  }

  # Compute statistics by time
  stats <- aggregate(
    as.formula(paste(value_column, "~ time")),
    data = ensemble_result,
    FUN = function(x) {
      c(
        n = length(x),
        mean = mean(x, na.rm = TRUE),
        sd = sd(x, na.rm = TRUE),
        min = min(x, na.rm = TRUE),
        max = max(x, na.rm = TRUE),
        quantile(x, probs = probs, na.rm = TRUE)
      )
    }
  )

  # Flatten the matrix columns
  stats_matrix <- stats[[value_column]]
  colnames(stats_matrix)[6:(5 + length(probs))] <- paste0("q", probs * 100)

  stats_df <- as.data.frame(stats_matrix)
  stats_df$time <- stats$time

  # Reorder columns
  time_col <- which(names(stats_df) == "time")
  stats_df <- stats_df[, c(time_col, setdiff(seq_len(ncol(stats_df)), time_col))]

  rownames(stats_df) <- NULL
  stats_df
}


#' @title Convert Ensemble to Wide Format
#'
#' @description Reshape ensemble data from long to wide format, with one column per model.
#' Useful for plotting and further analysis.
#'
#' @param ensemble_result Data frame returned by \code{\link{pv_dc_ensemble}} or
#'   \code{\link{pv_power_ensemble}}.
#' @param value_column Column name to reshape (default: "P_dc" for DC ensemble,
#'   "P_ac" for power ensemble). Set to NULL to keep all value columns.
#'
#' @return Data frame in wide format with one row per timestep and one column per model.
#'
#' @examples
#' \dontrun{
#' ens <- pv_dc_ensemble(...)
#' wide <- ensemble_wide(ens, "P_dc")
#' }
#'
#' @export
#'
ensemble_wide <- function(ensemble_result, value_column = NULL) {
  # Auto-detect value column if not specified
  if (is.null(value_column)) {
    # Keep all value columns in wide format
    value_columns <- c("P_dc", "P_ac", "G_poa", "T_cell", "iam")
    value_columns <- intersect(value_columns, names(ensemble_result))
  } else {
    value_columns <- value_column
  }

  # Keep only time, model, and value columns
  id_cols <- c("time", "model", "transposition", "decomposition", "cell_temp", "skoplaki_variant")
  keep_cols <- intersect(c("time", value_columns), names(ensemble_result))

  # Reshape each value column
  result <- ensemble_result[, c("time", "model", keep_cols), drop = FALSE]

  if (length(value_columns) == 1) {
    wide <- reshape(result[, c("time", "model", value_columns)],
                    idvar = "time",
                    timevar = "model",
                    direction = "wide")
    names(wide) <- gsub(paste0("^", value_columns, "\\."), "", names(wide))
  } else {
    # Use reshape2 or tidyr if available, otherwise base R
    if (requireNamespace("reshape2", quietly = TRUE)) {
      wide <- reshape2::dcast(result, time ~ model, value.var = value_columns)
    } else if (requireNamespace("tidyr", quietly = TRUE)) {
      wide_df <- result
      for (col in value_columns) {
        wide_df <- tidyr::pivot_wider(
          wide_df,
          names_from = "model",
          values_from = all_of(col),
          names_prefix = paste0(col, "_")
        )
      }
      wide <- wide_df
    } else {
      # Base R approach - do each column separately
      wide_list <- lapply(value_columns, function(col) {
        reshape(result[, c("time", "model", col), drop = FALSE],
               idvar = "time",
               timevar = "model",
               direction = "wide",
               v.names = col)
      })
      # Merge all
      wide <- wide_list[[1]]
      if (length(wide_list) > 1) {
        for (w in wide_list[-1]) {
          wide <- merge(wide, w, by = "time")
        }
      }
    }
  }

  rownames(wide) <- NULL
  wide
}


#' @title Rank Models by Performance
#'
#' @description Rank ensemble models by various performance metrics.
#'
#' @param ensemble_result Data frame returned by \code{\link{pv_dc_ensemble}} or
#'   \code{\link{pv_power_ensemble}}.
#' @param value_column Column name to analyze (default: "P_dc" for DC ensemble,
#'   "P_ac" for power ensemble).
#' @param metric Ranking metric: "mean" (average value), "max" (maximum value),
#'   "min" (minimum value), or "variance" (consistency).
#'
#' @return Data frame with model rankings.
#'
#' @examples
#' \dontrun{
#' ens <- pv_power_ensemble(...)
#' rankings <- ensemble_rank(ens, "P_ac", "mean")
#' }
#'
#' @export
#'
ensemble_rank <- function(ensemble_result, value_column = NULL, metric = c("mean", "max", "min", "variance")) {
  metric <- match.arg(metric)

  # Auto-detect value column if not specified
  if (is.null(value_column)) {
    if ("P_ac" %in% names(ensemble_result)) {
      value_column <- "P_ac"
    } else if ("P_dc" %in% names(ensemble_result)) {
      value_column <- "P_dc"
    } else {
      stop("Could not auto-detect value column. Please specify value_column.")
    }
  }

  if (!value_column %in% names(ensemble_result)) {
    stop(paste("Column", value_column, "not found in ensemble_result."))
  }

  # Compute metric by model
  model_stats <- aggregate(
    as.formula(paste(value_column, "~ model")),
    data = ensemble_result,
    FUN = function(x) {
      switch(metric,
             mean = mean(x, na.rm = TRUE),
             max = max(x, na.rm = TRUE),
             min = min(x, na.rm = TRUE),
             variance = var(x, na.rm = TRUE))
    }
  )

  names(model_stats)[[2]] <- metric

  # Add model metadata
  model_info <- unique(ensemble_result[, c("model", "transposition", "decomposition",
                                           "cell_temp", "skoplaki_variant")])
  model_stats <- merge(model_stats, model_info, by = "model")

  # Rank by metric (descending for mean/max, ascending for min/variance)
  if (metric %in% c("mean", "max")) {
    model_stats$rank <- rank(-model_stats[[metric]])
  } else {
    model_stats$rank <- rank(model_stats[[metric]])
  }

  # Sort by rank
  model_stats <- model_stats[order(model_stats$rank), ]
  rownames(model_stats) <- NULL

  model_stats[, c("rank", "model", metric, "transposition", "decomposition",
                  "cell_temp", "skoplaki_variant")]
}


#' @title Ensemble Spread Metrics
#'
#' @description Compute ensemble spread/uncertainty metrics for each timestep.
#'
#' @param ensemble_result Data frame returned by \code{\link{pv_dc_ensemble}} or
#'   \code{\link{pv_power_ensemble}}.
#' @param value_column Column name to analyze (default: "P_dc" for DC ensemble,
#'   "P_ac" for power ensemble).
#'
#' @return Data frame with spread metrics for each timestep.
#'
#' @examples
#' \dontrun{
#' ens <- pv_dc_ensemble(...)
#' spread <- pv_spread(ens, "P_dc")
#' }
#'
#' @export
#'
pv_spread <- function(ensemble_result, value_column = NULL) {
  # Auto-detect value column if not specified
  if (is.null(value_column)) {
    if ("P_ac" %in% names(ensemble_result)) {
      value_column <- "P_ac"
    } else if ("P_dc" %in% names(ensemble_result)) {
      value_column <- "P_dc"
    } else {
      stop("Could not auto-detect value column. Please specify value_column.")
    }
  }

  if (!value_column %in% names(ensemble_result)) {
    stop(paste("Column", value_column, "not found in ensemble_result."))
  }

  # Compute spread metrics by time
  spread <- aggregate(
    as.formula(paste(value_column, "~ time")),
    data = ensemble_result,
    FUN = function(x) {
      x <- x[!is.na(x)]
      if (length(x) < 2) {
        return(c(
          range = NA,
          iqr = NA,
          cv = NA,
          n_models = length(x)
        ))
      }
      c(
        range = max(x) - min(x),
        iqr = IQR(x),
        cv = ifelse(mean(x) != 0, sd(x) / abs(mean(x)), NA),
        n_models = length(x)
      )
    }
  )

  spread_matrix <- spread[[value_column]]
  spread_df <- as.data.frame(spread_matrix)
  spread_df$time <- spread$time

  # Reorder columns
  time_col <- which(names(spread_df) == "time")
  spread_df <- spread_df[, c(time_col, setdiff(seq_len(ncol(spread_df)), time_col))]

  rownames(spread_df) <- NULL
  names(spread_df) <- c("time", "range", "iqr", "cv", "n_models")

  spread_df
}


#' @title Prepare Ensemble Data for Plotting
#'
#' @description Transform ensemble data into a format suitable for common plotting libraries.
#'
#' @param ensemble_result Data frame returned by \code{\link{pv_dc_ensemble}} or
#'   \code{\link{pv_power_ensemble}}.
#' @param value_column Column name to plot (default: "P_dc" for DC ensemble,
#'   "P_ac" for power ensemble).
#' @param summarize Whether to return summary statistics instead of raw data (default FALSE).
#'
#' @return Data frame prepared for plotting with ggplot2 or base graphics.
#'
#' @examples
#' \dontrun{
#' ens <- pv_dc_ensemble(...)
#' plot_data <- ensemble_plot_data(ens, "P_dc")
#'
#' # With ggplot2
#' ggplot(plot_data, aes(x = time, y = P_dc, color = model)) +
#'   geom_line() +
#'   theme_minimal()
#' }
#'
#' @export
#'
ensemble_plot_data <- function(ensemble_result, value_column = NULL, summarize = FALSE) {
  # Auto-detect value column if not specified
  if (is.null(value_column)) {
    if ("P_ac" %in% names(ensemble_result)) {
      value_column <- "P_ac"
    } else if ("P_dc" %in% names(ensemble_result)) {
      value_column <- "P_dc"
    } else {
      stop("Could not auto-detect value column. Please specify value_column.")
    }
  }

  if (!value_column %in% names(ensemble_result)) {
    stop(paste("Column", value_column, "not found in ensemble_result."))
  }

  if (summarize) {
    # Return summary statistics for plotting
    summary_stats <- ensemble_summary(ensemble_result, value_column)
    return(summary_stats)
  }

  # Return raw data with necessary columns
  plot_cols <- c("time", "model", value_column)
  result <- ensemble_result[, intersect(plot_cols, names(ensemble_result)), drop = FALSE]

  result
}
