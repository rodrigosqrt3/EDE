#' Extinction date estimate
#'
#' Common S3 result class returned by every estimator in this package.
#'
#' Persistence tests ([solow1993()], [solow2005()], [mcinerny2006()],
#' [burgman1995()], [solow_roberts2003()] and [jaric2010()]) do not define a
#' point estimate. Inverting the test gives a one-sided `1 - alpha` upper
#' confidence bound for the extinction time, which is returned in `upper` with
#' `estimate = NA`. [solow1993b()] also returns the minimum variance unbiased
#' estimate of its model in `estimate`.
#'
#' @format A list with components:
#' \describe{
#'   \item{estimate}{Point estimate, or `NA` if not defined for this method.}
#'   \item{lower, upper}{Interval bounds, or `NA` if not defined. An upper
#'     bound that does not exist is returned as `Inf`.}
#'   \item{interval_type}{Either `"two-sided"` or `"one-sided"`.}
#'   \item{method}{Character string identifying the estimator.}
#'   \item{alpha}{Significance level used to compute the estimate.}
#' }
#' Bayesian estimators ([strauss1989_bayes()], [weiss_marshall1999()]) return
#' posterior probability (credible) intervals and set `interval_kind` to
#' `"credible"`. Some estimators append method-specific components,
#' documented on their help pages.
#' @name ede_estimate
NULL

#' @keywords internal
#' @noRd
new_ede_estimate <- function(estimate, lower = NA_real_, upper = NA_real_,
                             method, alpha, interval_type = "two-sided") {
  structure(
    list(
      estimate = estimate,
      lower = lower,
      upper = upper,
      method = method,
      alpha = alpha,
      interval_type = interval_type
    ),
    class = "ede_estimate"
  )
}

#' Result of inverting a persistence test
#'
#' `bound` is the earliest candidate extinction time at which persistence is
#' rejected at level `alpha`, i.e. a one-sided `1 - alpha` upper confidence
#' bound. Bounds that do not exist or lie beyond `test_year` give `NA` with a
#' warning.
#' @keywords internal
#' @noRd
test_bound_result <- function(bound, test_year, method, alpha) {
  if (is.na(bound) || bound > test_year) {
    warning(
      "p-value never falls to alpha before `test_year`; ",
      "returning NA.", call. = FALSE
    )
    bound <- NA_real_
  }
  new_ede_estimate(NA_real_, upper = bound, method = method, alpha = alpha,
                   interval_type = "one-sided")
}

#' Validate `test_year`, and if given, compare it with the last sighting
#' @keywords internal
#' @noRd
check_test_year <- function(test_year, last_sight = NULL) {
  if (missing(test_year) || !is.numeric(test_year) ||
      length(test_year) != 1L || is.na(test_year)) {
    stop("`test_year` must be supplied as a number.", call. = FALSE)
  }
  if (!is.null(last_sight) && test_year <= last_sight) {
    stop("`test_year` must be later than the last sighting.", call. = FALSE)
  }
  invisible(test_year)
}

#' Integer-step candidate times after the last sighting, up to `test_year`
#' @keywords internal
#' @noRd
candidate_times <- function(last_sight, test_year) {
  if (test_year < last_sight + 1) return(numeric(0))
  seq(last_sight + 1, test_year)
}

#' @export
print.ede_estimate <- function(x, ...) {
  cat(sprintf("<%s>\n", x$method))
  # persistence tests define only an upper bound, possibly NA
  bound_only <- is.na(x$estimate) && is.na(x$lower) &&
    identical(x$interval_type, "one-sided")
  if (!bound_only) {
    cat(sprintf("  estimate: %s\n", format(round(x$estimate, 2))))
  }
  credible <- identical(x$interval_kind, "credible")
  one_sided <- identical(x$interval_type, "one-sided")
  if (one_sided && is.na(x$lower)) {
    cat(sprintf(
      "  %g%% upper %s bound: %s\n",
      100 * (1 - x$alpha), if (credible) "credible" else "confidence",
      format(round(x$upper, 2))
    ))
    return(invisible(x))
  }
  if (!is.na(x$lower) || !is.na(x$upper)) {
    interval_label <- paste0(if (one_sided) "one-sided " else "",
                             if (credible) "credible interval" else "CI")
    cat(sprintf(
      "  %g%% %s: [%s, %s]\n",
      100 * (1 - x$alpha), interval_label,
      format(round(x$lower, 2)), format(round(x$upper, 2))
    ))
  }
  invisible(x)
}
