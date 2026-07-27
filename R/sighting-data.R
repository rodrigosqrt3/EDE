#' Construct a validated sighting record
#'
#' Builds the common input object used by every estimator in the package:
#' a time-ordered table of sighting counts, checked for the conditions each
#' estimator in the sighting-record literature assumes (numeric, non-negative,
#' no duplicated times).
#'
#' @param data A data frame or matrix. By default the first column is read
#'   as time (e.g. year) and the second as the number of sightings recorded
#'   at that time.
#' @param time_col,count_col Column name or position for time and sighting
#'   count.
#'
#' @return An object of class `sighting_data`: a data frame with columns
#'   `time` and `count`, sorted by time.
#'
#' @export
sighting_data <- function(data, time_col = 1L, count_col = 2L) {
  if (missing(data)) stop("`data` is required.", call. = FALSE)
  if (!is.data.frame(data) && !is.matrix(data)) {
    stop("`data` must be a data.frame or matrix.", call. = FALSE)
  }
  data <- as.data.frame(data)
  if (ncol(data) < 2L) stop("`data` needs at least two columns.", call. = FALSE)

  time <- data[[time_col]]
  count <- data[[count_col]]

  if (!is.numeric(time) || !is.numeric(count)) {
    stop("time and count columns must be numeric.", call. = FALSE)
  }
  if (anyNA(time) || anyNA(count)) {
    stop("time and count columns cannot contain NA.", call. = FALSE)
  }
  if (any(time < 0) || any(count < 0)) {
    stop("time and count cannot be negative.", call. = FALSE)
  }
  if (anyDuplicated(time)) {
    stop(
      "`time` has duplicated values. Aggregate sightings per time unit ",
      "before calling sighting_data().",
      call. = FALSE
    )
  }
  if (max(count) > max(time)) {
    warning(
      "sighting counts are larger than any time value; check that ",
      "time_col/count_col point at the right columns.",
      call. = FALSE
    )
  }

  ord <- order(time)
  out <- data.frame(time = time[ord], count = count[ord])
  class(out) <- c("sighting_data", class(out))
  out
}

#' Expand a sighting record to one row per time unit up to `end_time`
#'
#' Internal helper used by likelihood-based estimators (Solow, Burgman) that
#' need the full sequence of zero/nonzero sighting counts between the first
#' sighting and a candidate test year, not just the nonzero rows.
#'
#' @param sd A [sighting_data] object.
#' @param end_time Last time unit to include.
#' @keywords internal
#' @noRd
expand_record <- function(sd, end_time) {
  sd <- sd[sd$count > 0, , drop = FALSE]
  if (nrow(sd) == 0L) stop("no sightings with count > 0.", call. = FALSE)
  if (end_time < sd$time[1]) {
    stop("`end_time` cannot precede the first sighting.", call. = FALSE)
  }
  if (end_time < sd$time[nrow(sd)]) {
    stop("`end_time`/`test_year` must be later than the last sighting.", call. = FALSE)
  }
  time <- seq(sd$time[1], end_time)
  count <- rep(0, length(time))
  count[match(sd$time, time)] <- sd$count
  data.frame(time = time, count = count)
}
