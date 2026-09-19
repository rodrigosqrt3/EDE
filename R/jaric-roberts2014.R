#' Reliability-adjusted extinction inference
#'
#' Modifies the constant-rate Solow model by assigning a reliability in
#' `[0, 1]` to each occupied sighting time. Implements equations 4, 6, and
#' 8-10 of Jaric & Roberts (2014), returning the reliability-adjusted point
#' estimate and upper confidence bound.
#'
#' @param sd A [sighting_data] object. This method currently requires binary
#'   counts (zero or one) because reliability is assigned to individual
#'   observations.
#' @param reliability Numeric reliabilities in `[0, 1]`, supplied either for
#'   every row of `sd` or for the positive-count rows only.
#' @param alpha Significance level, in (0, 1), for the upper confidence bound.
#'
#' @return An [ede_estimate] object.
#'
#' @references
#' Jaric, I., & Roberts, D. L. (2014). Accounting for observation reliability
#'   when inferring extinction based on sighting records. Biodiversity and
#'   Conservation, 23(11), 2801-2815.
#'
#' @export
jaric_roberts2014 <- function(sd, reliability, alpha = 0.05) {
  stopifnot(inherits(sd, "sighting_data"))
  if (alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0, 1).", call. = FALSE)
  if (missing(reliability) || !is.numeric(reliability) || anyNA(reliability) ||
      any(!is.finite(reliability)) || any(reliability < 0 | reliability > 1)) {
    stop("`reliability` must contain finite values in [0, 1].", call. = FALSE)
  }
  if (any(sd$count > 1)) {
    stop("Jaric & Roberts (2014) currently requires binary sighting counts.",
         call. = FALSE)
  }

  occupied <- which(sd$count > 0)
  if (length(reliability) == nrow(sd)) {
    event_reliability <- reliability[occupied]
  } else if (length(reliability) == length(occupied)) {
    event_reliability <- reliability
  } else {
    stop("`reliability` must have one value per row or per occupied sighting time.",
         call. = FALSE)
  }

  event_times <- sd$time[occupied]
  origin <- min(sd$time)
  used <- event_times > origin
  event_times <- event_times[used] - origin
  event_reliability <- event_reliability[used]
  if (length(event_times) == 0L) {
    stop("at least one sighting is required after the observation origin.",
         call. = FALSE)
  }

  r <- sum(event_reliability)
  if (r <= 0) {
    stop("the sum of reliabilities after the observation origin must be positive.",
         call. = FALSE)
  }

  n <- length(event_times)
  last_valid_probability <- vapply(seq_len(n), function(i) {
    later_false <- if (i == n) 1 else prod(1 - event_reliability[(i + 1L):n])
    event_reliability[i] * later_false
  }, numeric(1))
  tr <- sum(event_times * last_valid_probability)

  estimate <- origin + ((r + 1) / r) * tr
  upper <- origin + tr / alpha^(1 / r)
  result <- new_ede_estimate(
    estimate,
    lower = NA_real_,
    upper = upper,
    method = "Jaric & Roberts (2014)",
    alpha = alpha,
    interval_type = "one-sided"
  )
  result$effective_sightings <- r
  result$effective_endpoint <- origin + tr
  result
}
