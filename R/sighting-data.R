#' Construct a validated sighting record
#'
#' Builds the common input object used by every estimator in the package:
#' a time-ordered table of sighting counts, checked for the conditions each
#' estimator in the sighting-record literature assumes (numeric, non-negative,
#' no duplicated times).
#'
#' @param data A data frame or matrix. By default the first column is read
#'   as time (e.g. year) and the second as the number of sightings recorded
#'   at that time. The earliest supplied time defines the observation origin;
#'   include an initial zero-count row when observation began before the first
#'   sighting.
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
  if (anyNA(time) || anyNA(count) || any(!is.finite(time)) || any(!is.finite(count))) {
    stop("time and count columns must contain finite, non-missing values.", call. = FALSE)
  }
  if (any(count < 0)) {
    stop("sighting counts cannot be negative.", call. = FALSE)
  }
  if (any(count != floor(count))) {
    stop("sighting counts must be whole numbers.", call. = FALSE)
  }
  if (anyDuplicated(time)) {
    stop(
      "`time` has duplicated values. Aggregate sightings per time unit ",
      "before calling sighting_data().",
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
  if (!any(sd$count > 0)) stop("no sightings with count > 0.", call. = FALSE)
  start_time <- min(sd$time)
  last_sighting <- max(sd$time[sd$count > 0])
  if (end_time < start_time) {
    stop("`end_time` cannot precede the first sighting.", call. = FALSE)
  }
  if (end_time < last_sighting) {
    stop("`end_time`/`test_year` must be later than the last sighting.", call. = FALSE)
  }
  time <- seq(start_time, end_time)
  count <- rep(0, length(time))
  count[match(sd$time, time)] <- sd$count
  data.frame(time = time, count = count)
}
