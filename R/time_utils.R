#' @title Time Zone Handling Utilities
#'
#' @description Internal functions for consistent timezone handling throughout
#' the pvflux package. All solar position calculations are performed in UTC
#' to ensure consistency, with timestamps converted back to their original
#' timezone in the output.
#'
#' @details
#' The pvflux package follows these timezone conventions:
#' \itemize{
#'   \item Input timestamps are converted to POSIXct if not already
#'   \item If no timezone is specified, UTC is assumed
#'   \item All internal calculations use UTC time
#'   \item Output timestamps are returned in the original input timezone
#' }
#'
#' This approach ensures that solar position calculations are always performed

#' consistently regardless of the input timezone, while preserving the user's
#' preferred time representation in the output.
#'
#' @name time_utils
#' @keywords internal
NULL


#' Prepare Time for UTC-based Calculations
#'
#' Converts input time to POSIXct in UTC for internal calculations.
#' Stores the original timezone for later restoration.
#'
#' @param time A vector of timestamps. Can be POSIXct, POSIXlt, character,
#'   or numeric (interpreted as seconds since epoch in UTC).
#'
#' @return A list with two elements:
#' \describe{
#'   \item{time_utc}{POSIXct vector in UTC timezone}
#'   \item{original_tz}{Character string of original timezone (or "UTC" if none)}
#' }
#'
#' @details
#' The function handles several input types:
#' \itemize{
#'   \item \strong{POSIXct with timezone}: Converts to UTC, preserves original tz
#'   \item \strong{POSIXct without timezone}: Assumes UTC
#'   \item \strong{POSIXlt}: Converts to POSIXct, then handles as above
#'   \item \strong{Character}: Parses with lubridate, assumes UTC if no tz info
#'   \item \strong{Numeric}: Interprets as seconds since 1970-01-01 00:00:00 UTC
#' }
#'
#' @examples
#' \dontrun{
#' # POSIXct with explicit timezone - converts to UTC
#' t1 <- as.POSIXct("2026-01-15 12:00", tz = "Africa/Johannesburg")
#' prep1 <- prepare_time_utc(t1)
#' # prep1$time_utc is 10:00 UTC (SAST is UTC+2)
#' # prep1$original_tz is "Africa/Johannesburg"
#'
#' # POSIXct with UTC - no conversion needed
#' t2 <- as.POSIXct("2026-01-15 12:00", tz = "UTC")
#' prep2 <- prepare_time_utc(t2)
#' # prep2$time_utc is 12:00 UTC
#' # prep2$original_tz is "UTC"
#'
#' # POSIXct without timezone - assumes UTC
#' t3 <- as.POSIXct("2026-01-15 12:00")
#' prep3 <- prepare_time_utc(t3)
#' # prep3$time_utc is 12:00 UTC
#' # prep3$original_tz is "UTC"
#' }
#'
#' @keywords internal
prepare_time_utc <- function(time) {
  # Handle POSIXlt by converting to POSIXct first

if (inherits(time, "POSIXlt")) {
    time <- as.POSIXct(time)
  }

  # Handle character input
  if (is.character(time)) {
    time <- lubridate::ymd_hms(time, quiet = TRUE)
    if (any(is.na(time))) {
      # Try other formats
      time_orig <- time
      time <- lubridate::parse_date_time(time_orig, orders = c("ymd HMS", "ymd HM", "ymd"))
    }
    # If still no tz, lubridate defaults to UTC
  }

  # Handle numeric input (seconds since epoch)
  if (is.numeric(time)) {
    time <- as.POSIXct(time, origin = "1970-01-01", tz = "UTC")
  }

  # Now time should be POSIXct
  if (!inherits(time, "POSIXct")) {
    stop("Could not convert 'time' to POSIXct. Please provide POSIXct, ",
         "POSIXlt, character, or numeric timestamps.")
  }

  # Get the original timezone
  original_tz <- attr(time, "tzone")

  # Handle empty or NULL timezone (system default) - treat as UTC
  if (is.null(original_tz) || length(original_tz) == 0 || original_tz[1] == "") {
    original_tz <- "UTC"
    # Set UTC explicitly so conversion works correctly
    time <- lubridate::force_tz(time, tzone = "UTC")
  } else {
    original_tz <- original_tz[1]  # Take first element if vector
  }

  # Convert to UTC for calculations
  if (original_tz != "UTC") {
    time_utc <- lubridate::with_tz(time, tzone = "UTC")
  } else {
    time_utc <- time
  }

  list(
    time_utc = time_utc,
    original_tz = original_tz
  )
}


#' Restore Time to Original Timezone
#'
#' Converts UTC timestamps back to the original timezone for output.
#'
#' @param time_utc POSIXct vector in UTC timezone.
#' @param original_tz Character string of the target timezone.
#'
#' @return POSIXct vector in the original timezone.
#'
#' @examples
#' \dontrun{
#' # After calculations in UTC, restore to original SAST timezone
#' time_utc <- as.POSIXct("2026-01-15 10:00", tz = "UTC")
#' time_sast <- restore_time_tz(time_utc, "Africa/Johannesburg")
#' # time_sast shows 12:00 SAST
#' }
#'
#' @keywords internal
restore_time_tz <- function(time_utc, original_tz) {
  if (original_tz == "UTC") {
    return(time_utc)
  }
  lubridate::with_tz(time_utc, tzone = original_tz)
}
