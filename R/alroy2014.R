#' Alroy (2014) "creeping-shadow-of-a-doubt" posterior probability of
#' extinction
#'
#' Posterior probability that a species is extinct, interval by interval after
#' its last sighting, by the iterative Bayesian method of Alroy (2014). The
#' record is a series of time units (years, say) in each of which the species
#' was either sighted (a presence) or not (an absence). The prior probability
#' of going extinct in one unit, `E`, comes from the length of the observed
#' range, by assuming that extinction is an exponential decay process with a
#' 50% chance that the species has gone extinct by the end of the range; the
#' sampling probability `ps` is the frequency of sightings within the range.
#' Bayes' theorem then gives the posterior `e[1]` for the first unit after the
#' last sighting, which is carried to the next unit as the prior, where
#' extinction might instead have happened:
#' `e[t] = (e[t-1] + (1 - e[t-1]) E) / ((e[t-1] + (1 - e[t-1]) E) +
#' (1 - e[t-1]) (1 - E) (1 - ps))` (equations 3 and 4 of the paper), starting
#' from `e[0] = 0` at the last sighting. The values are Alroy's "shadow values".
#'
#' The range `r` is the number of time units from `base` (by default the first
#' sighting) to the last sighting, both included, and the sighting
#' probability uses the `np` units with sightings and the `na` units without
#' them between the first and the last sighting, `ps = (np - 2) / (np - 2 +
#' na)` (the first and last presences are required to exist). With
#' `adjust_range = TRUE` the range is replaced by `r (np + 1) / (np - 1)`, the
#' correction of Strauss & Sadler (1989) used in the simulations of the paper.
#' The paper's remarks on how to choose the range apply: the species may be
#' older than the record (`truncated = TRUE` uses `ps = (np - 1) / (np - 1 +
#' na)`), and for a recent extinction the range should start at the first date
#' at which extinction was plausible rather than at the first sighting, which
#' `base` allows; results sensitive to that choice should be treated with
#' caution. With `run_based = TRUE` the sighting probability is conditioned on
#' the state of the preceding unit, `ps` given a presence for the first unit and
#' given an absence afterwards (the "run-based" version of the paper, which it
#' recommends when the sample is large enough).
#'
#' Following the paper, the best estimate of the extinction time is the first
#' unit at which the shadow value reaches 0.5, and credible limits are the
#' first units at which it exceeds `1 - alpha` (one-sided) or `alpha / 2` and
#' `1 - alpha / 2` (two-sided). A species with fewer than three units with
#' sightings cannot be analysed. The versions of the paper that use sampling
#' controls (several sampling events per unit) and the maximum-likelihood
#' alternative of its Alternative Method section are not implemented.
#'
#' @param sd A [sighting_data] object whose times are whole numbers (time
#'   units); a unit is a presence if its count is positive.
#' @param test_year Latest time unit for which a shadow value is computed, a
#'   whole number later than the last sighting.
#' @param alpha Significance level, in (0, 1), for the credible limits.
#' @param interval `"one-sided"` for an upper credible limit, or `"two-sided"`
#'   for an equal-tailed credible interval.
#' @param base First time unit at which extinction was plausible; defaults to
#'   the first sighting and cannot be later than it.
#' @param adjust_range Correct the range for sighting failures
#'   (`r (np + 1) / (np - 1)`).
#' @param truncated If `TRUE`, the first sighting is not the first appearance of
#'   the species, which changes the sampling probability (see above).
#' @param run_based Use the sampling probabilities conditioned on the preceding
#'   unit.
#' @param data_out If `TRUE`, return the shadow values as a data frame instead
#'   of the estimate.
#'
#' @return An [ede_estimate] object with `estimate` the first time unit at which
#'   the shadow value is at least 0.5 and the credible limits in `lower` and
#'   `upper` (`NA`, with a warning, when not reached by `test_year`), and the
#'   components `shadow` (a data frame with `time` and `e`), `E`, `ps`, `r`,
#'   `n_presence` and `n_absence`. With `data_out = TRUE` the data frame of
#'   shadow values.
#'
#' @references
#' Alroy, J. (2014). A simple Bayesian method of inferring extinction.
#'   Paleobiology, 40(4), 584-607. doi:10.1666/13074
#'
#' Strauss, D., & Sadler, P. M. (1989). Classical confidence intervals and
#'   Bayesian probability estimates for ends of local taxon ranges.
#'   Mathematical Geology, 21(4), 411-427.
#'
#' @seealso [weiss_marshall1999()] and [strauss1989_bayes()] for other Bayesian
#'   estimates of the end of a range.
#'
#' @examples
#' # a species seen in 14 of 50 units (about once every four), then not seen
#' times <- c(0, 4, 7, 11, 15, 18, 22, 26, 29, 33, 37, 41, 45, 49)
#' sd <- sighting_data(data.frame(time = times, count = 1))
#' alroy2014(sd, test_year = 80)
#' head(alroy2014(sd, test_year = 80, data_out = TRUE))
#' alroy2014(sd, test_year = 80, run_based = TRUE, interval = "two-sided")
#'
#' @export
alroy2014 <- function(sd, test_year, alpha = 0.05,
                      interval = c("one-sided", "two-sided"), base = NULL,
                      adjust_range = TRUE, truncated = FALSE,
                      run_based = FALSE, data_out = FALSE) {
  stopifnot(inherits(sd, "sighting_data"))
  check_alpha(alpha)
  interval <- match.arg(interval)
  pos <- sd$count > 0
  if (!any(pos)) stop("no sightings with count > 0.", call. = FALSE)
  tp <- sd$time[pos]
  if (any(tp != round(tp))) {
    stop("`sd$time` must be whole numbers (time units).", call. = FALSE)
  }
  first <- min(tp)
  last <- max(tp)
  check_test_year(test_year, last)
  if (test_year != round(test_year)) {
    stop("`test_year` must be a whole number.", call. = FALSE)
  }
  present <- unique(tp)
  np <- length(present)
  if (np < 3L) {
    stop("Alroy (2014) needs at least 3 time units with sightings.",
         call. = FALSE)
  }
  if (is.null(base)) base <- first
  if (!is.numeric(base) || length(base) != 1L || is.na(base) ||
      base != round(base) || base > first) {
    stop("`base` must be a whole number not later than the first sighting.",
         call. = FALSE)
  }

  grid <- first:last
  x <- grid %in% present
  na <- sum(!x)
  r <- last - base + 1
  r_adj <- if (isTRUE(adjust_range)) r * (np + 1) / (np - 1) else r
  big_e <- 1 - exp(log(0.5) / r_adj) # 1 - exp(-lambda), lambda = -log(.5) / r

  ps <- if (isTRUE(truncated)) (np - 1) / (np - 1 + na) else
    (np - 2) / (np - 2 + na)
  ps_first <- ps_next <- ps
  if (isTRUE(run_based)) {
    prev <- x[-length(x)]
    cur <- x[-1L]
    npp <- sum(prev & cur)
    npa <- sum(prev & !cur)
    nap <- sum(!prev & cur)
    naa <- sum(!prev & !cur)
    ps_first <- if (npa == 0) npp / (npp + 1) else npp / (npp + npa)
    ps_next <- if (naa == 0) nap / (nap + 1) else nap / (nap + naa)
  }

  n <- test_year - last
  e <- numeric(n)
  prior <- 0
  for (t in seq_len(n)) {
    s_t <- if (t == 1L) ps_first else ps_next
    a <- prior + (1 - prior) * big_e
    e[t] <- a / (a + (1 - prior) * (1 - big_e) * (1 - s_t))
    prior <- e[t]
  }
  shadow <- data.frame(time = last + seq_len(n), e = e)
  if (isTRUE(data_out)) return(shadow)

  reach <- function(p) {
    hit <- which(e >= p)
    if (length(hit)) shadow$time[hit[1L]] else NA_real_
  }
  est <- reach(0.5)
  if (interval == "one-sided") {
    lower <- NA_real_
    upper <- reach(1 - alpha)
  } else {
    lower <- reach(alpha / 2)
    upper <- reach(1 - alpha / 2)
  }
  if (is.na(est) || is.na(upper)) {
    warning("the shadow value does not reach the required level before ",
            "`test_year`; the corresponding limit is NA.", call. = FALSE)
  }
  out <- new_ede_estimate(est, lower, upper, method = "Alroy (2014)",
                          alpha = alpha, interval_type = interval)
  out$interval_kind <- "credible"
  out$shadow <- shadow
  out$E <- big_e
  out$ps <- if (isTRUE(run_based)) c(first = ps_first, next_unit = ps_next) else ps
  out$r <- r
  out$n_presence <- np
  out$n_absence <- na
  out
}
