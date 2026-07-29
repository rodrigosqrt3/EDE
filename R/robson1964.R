#' Robson & Whitlock (1964) truncation point estimator
#'
#' Extrapolates the extinction date from the gap between the two most recent
#' sighting times, under the assumption that the sighting process near the
#' true endpoint behaves like the tail of a uniform record process.
#'
#' @param sd A [sighting_data] object.
#' @param alpha Significance level, in (0, 1). Not a coverage-calibrated CI
#'   here: it directly scales the extrapolated gap, following the original
#'   formula, so there is no `lower`/`upper` in the returned estimate.
#'
#' @return An [ede_estimate] object.
#'
#' @references
#' Robson, D. S., & Whitlock, J. H. (1964). Estimation of a truncation
#'   point. Biometrika, 51(1/2), 33-39.
#'
#' @export
robson1964 <- function(sd, alpha = 0.05) {
  stopifnot(inherits(sd, "sighting_data"))
  if (alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0, 1).", call. = FALSE)

  times <- sort(unique(sd$time[sd$count > 0]))
  n <- length(times)
  if (n < 2) {
    stop("Robson & Whitlock (1964) needs at least 2 distinct sighting times with count > 0.", call. = FALSE)
  }

  gap <- times[n] - times[n - 1]
  estimate <- times[n] + gap * (1 - alpha) / alpha

  new_ede_estimate(estimate, lower = NA_real_, upper = NA_real_, method = "Robson & Whitlock (1964)", alpha = alpha)
}
