#' Marshall (1997) generalized confidence interval for a range end-point
#'
#' Confidence interval on the end-point of a stratigraphic (or sighting) range
#' when the probability of recovering a find varies along the axis according
#' to a user-supplied fossil recovery potential `f(t)`. With `H` distinct find
#' horizons, the observed range from the first (FAD) to the last (LAD) find
#' carries recovery mass `R* = integral of f from FAD to LAD` (equation 10 of
#' Marshall, 2010). The `C = 1 - alpha` interval extends beyond the LAD by the
#' distance `r` whose recovery mass `r* = integral of f from LAD to LAD + r`
#' (equation 11) equals `R* [(1 - C)^(-1 / (H - 1)) - 1]` (equation 9, the
#' generalization of the classical interval, equation 4). If even the whole
#' recovery mass beyond the LAD (up to `max_time`) is smaller than that, the
#' interval is unbounded and `upper = Inf`, as in the short half-life case of
#' Marshall (2010, Case 3). With constant `f` the method reduces to the
#' classical interval of [strauss1989()].
#'
#' The confidence that the end-point lies before a candidate time `T`,
#' `C = 1 - [R* / (R* + r*)]^(H - 1)` with `r*` the recovery mass from the LAD
#' to `T` (equations 3 and 6 with `R*`, `r*` in place of the range and its
#' extension), is returned as `confidence` when `test_year` is given.
#'
#' Each distinct time with a positive count is one find horizon; counts are
#' ignored and zero-count rows are not used, as in [strauss1989()]. The end
#' point estimated is always the one at the largest `time`: to place an
#' interval on a last appearance dated in Ma, code `time = -age`; for a first
#' appearance, code `time = -height` (or `time = age`). `recovery` must be
#' written on the same coded axis.
#'
#' The formulas are from Marshall (2010, Box 3), who attributes the method to
#' Marshall (1997); they have been checked against the worked example of
#' Marshall (2010, section 4.1), not against Marshall (1997) itself. The
#' integrals are computed numerically with [stats::integrate()], so a
#' discontinuous (e.g. stepwise) recovery potential should be supplied with
#' breakpoints smoothed or with care.
#'
#' @param sd A [sighting_data] object.
#' @param recovery A vectorized function of time giving the (relative) fossil
#'   recovery potential, non-negative. Only ratios of its integrals matter, so
#'   it need not be normalized. The default (constant) gives the classical
#'   interval.
#' @param alpha Significance level, in (0, 1); the confidence level is
#'   `1 - alpha`.
#' @param test_year Optional time after the last find at which to compute the
#'   confidence that the end-point had already been reached.
#' @param max_time Upper limit of the axis on which `recovery` is defined
#'   (e.g. the top of the section). Defaults to `Inf`.
#'
#' @return An [ede_estimate] object with `estimate = NA`, `lower = NA`,
#'   `upper` = the end of the one-sided `1 - alpha` interval (possibly `Inf`)
#'   and `interval_type = "one-sided"`, plus components `last_sighting` (the
#'   LAD), `recovery_range` (`R*`), `recovery_extension` (the `r*` required)
#'   and, if `test_year` is given, `confidence`.
#'
#' @references
#' Marshall, C. R. (1997). Confidence intervals on stratigraphic ranges with
#'   non-random distributions of fossil horizons. Paleobiology, 23(2),
#'   165-173.
#'
#' Marshall, C. R. (2010). Using confidence intervals to quantify the
#'   uncertainty in the end-points of stratigraphic ranges. In J. Alroy & G.
#'   Hunt (Eds.), Quantitative Methods in Paleobiology (pp. 291-316). The
#'   Paleontological Society Papers, 16.
#'
#' @export
marshall1997 <- function(sd, recovery = function(t) rep(1, length(t)),
                         alpha = 0.05, test_year = NULL, max_time = Inf) {
  stopifnot(inherits(sd, "sighting_data"))
  if (alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0, 1).", call. = FALSE)
  if (!is.function(recovery)) stop("`recovery` must be a function of time.", call. = FALSE)

  times <- sort(unique(sd$time[sd$count > 0]))
  h <- length(times)
  if (h < 2) {
    stop("Marshall (1997) needs at least 2 distinct sighting times with count > 0.",
         call. = FALSE)
  }
  fad <- times[1L]
  lad <- times[h]
  if (!is.numeric(max_time) || length(max_time) != 1L || is.na(max_time) ||
      max_time <= lad) {
    stop("`max_time` must be a number later than the last sighting.", call. = FALSE)
  }
  if (!is.null(test_year)) check_test_year(test_year, lad)

  f <- function(t) {
    v <- recovery(t)
    if (length(v) != length(t) || !is.numeric(v) || anyNA(v) || any(v < 0)) {
      stop("`recovery` must return one non-negative number per time.", call. = FALSE)
    }
    v
  }
  mass <- function(from, to) {
    if (to <= from) return(0)
    stats::integrate(f, from, to, subdivisions = 1000L,
                     rel.tol = 1e-10, stop.on.error = TRUE)$value
  }

  r_star <- mass(fad, lad)
  if (r_star <= 0) {
    stop("the recovery potential integrates to zero over the observed range.",
         call. = FALSE)
  }
  ext <- r_star * (alpha^(-1 / (h - 1)) - 1)

  # a divergent integral (e.g. constant f on an unbounded axis) means the
  # tail mass is unlimited
  tail_mass <- tryCatch(mass(lad, max_time), error = function(e) {
    if (is.infinite(max_time) && grepl("divergent", conditionMessage(e))) Inf
    else stop(e)
  })
  if (tail_mass < ext) {
    upper <- Inf
  } else {
    # bracket the root of mass(lad, lad + r) = ext
    step <- (lad - fad)
    hi <- min(lad + step, max_time)
    while (mass(lad, hi) < ext && hi < max_time) {
      step <- step * 2
      hi <- min(lad + step, max_time)
    }
    upper <- stats::uniroot(function(x) mass(lad, x) - ext,
                            c(lad, hi), tol = 1e-10 * max(1, abs(hi)))$root
  }

  result <- new_ede_estimate(NA_real_, upper = upper,
                             method = "Marshall (1997)", alpha = alpha,
                             interval_type = "one-sided")
  result$last_sighting <- lad
  result$recovery_range <- r_star
  result$recovery_extension <- ext
  if (!is.null(test_year)) {
    rt <- mass(lad, min(test_year, max_time))
    result$confidence <- 1 - (r_star / (r_star + rt))^(h - 1)
  }
  result
}
