#' Rout et al. (2009) Bayes factor with a declining sighting rate
#'
#' Bayesian test of persistence (e.g. of an invasive species after an
#' eradication programme) in which the pre-extinction sighting rate declines
#' as `lambda(s) = m * s^(-a)` (equation 6 of Rout, Salomon & McCarthy,
#' 2009), with `m` a constant and `0 <= a < 1`. Integrating over `m` (with the
#' prior of Solow, 1993a) gives the Bayes factor in favour of persistence to
#' `S`:
#' `B = (n(1 - a) - 1) / ((S / s[n])^(n(1 - a) - 1) - 1)`
#' (Rout et al., 2009, unnumbered equation after equation 6; derivation in
#' their Appendix S1). With `a = 0` it is the Bayes factor of Solow (1993a)
#' (Rout et al., equation 2), and this function then equals
#' [solow1993_bayes()]. When `n(1 - a) = 1` the limit `1 / log(S / s[n])` is
#' used. The posterior probability of persistence is
#' `1 / (1 + (1 - prior) / (prior * B))` (equation 1).
#'
#' Rout et al. (2009) measure time in surveys: `S` surveys have been made and
#' the species was seen in `n` of them, the last being survey `s[n]`. Here
#' times are measured from the observation origin, the earliest time in `sd`,
#' and sighting events are counted as in [solow1993()]: a sighting at the
#' origin defines it and is not counted. To follow the paper, supply one row
#' per survey with binary counts and a leading zero-count row for survey 0.
#'
#' The paper's model has the rate declining with the number of surveys, so
#' the result depends on where the origin is placed; this is a property of
#' the model, not of the implementation.
#'
#' @inheritParams solow1993_bayes
#' @param a Rate of decline of the sighting rate, in `[0, 1)`. Rout et al.
#'   (2009) estimated `a = 0.1` for *Helenium amarum*.
#' @param test_year Candidate time `S` (e.g. the current survey number). Must
#'   be later than the last sighting.
#'
#' @return A list of class `ede_bayes` with components `bayes_factor`,
#'   `posterior`, `prior`, `test_year`, `method` and `a`; or (if
#'   `data_out = TRUE`) a data frame with columns `time`, `bayes_factor` and
#'   `posterior`.
#'
#' @seealso [rout2009_absences()] for the number of absent surveys after
#'   which eradication should be declared.
#'
#' @references
#' Rout, T. M., Salomon, Y., & McCarthy, M. A. (2009). Using sighting
#'   records to declare eradication of an invasive species. Journal of
#'   Applied Ecology, 46(1), 110-117.
#'
#' Solow, A. R. (1993a). Inferring extinction from sighting data. Ecology,
#'   74(3), 962-964.
#'
#' @export
rout2009 <- function(sd, test_year, a = 0, prior = 0.5, data_out = FALSE) {
  info <- rout2009_setup(sd, a, prior)
  check_test_year(test_year, info$last_sight)
  posterior <- function(b) 1 / (1 + (1 - prior) / (prior * b))
  bf <- function(t) rout2009_bf(info$n, info$tn, t - info$origin, a)

  if (data_out) {
    times <- candidate_times(info$last_sight, test_year)
    b <- vapply(times, bf, numeric(1))
    return(data.frame(time = times, bayes_factor = b, posterior = posterior(b)))
  }

  b <- bf(test_year)
  structure(
    list(bayes_factor = b, posterior = posterior(b), prior = prior,
         test_year = test_year,
         method = sprintf("Rout et al. (2009) Bayes factor (a = %g)", a),
         a = a),
    class = "ede_bayes"
  )
}

#' Optimal number of absent surveys before declaring eradication
#'
#' Decision rule of Rout, Salomon & McCarthy (2009): after `d` consecutive
#' surveys without a sighting, the net expected cost of stopping, in units
#' of the cost of escape and damage `C[e]`, is
#' `NEC(d) / C[e] = (d - 1) R + p(extant | d)` (equation 3, the "rule of
#' thumb"), where `R = C[s] / C[e]` is the ratio of the cost of one survey to
#' the cost of escape, and `p(extant | d)` is the posterior probability of
#' persistence of [rout2009()] with `S = s[n] + d`. The optimal number of
#' absent surveys `d*` minimises `NEC(d)` over `d = 0, 1, ..., max_absent`.
#'
#' For the constant-rate model (`a = 0`) the closed-form approximation
#' `d* ~ s[n] ((n - 1)^2 prior / (s[n] (1 - prior) R))^(1 / n) - s[n]`
#' (equation 5) is also returned. The paper gives no approximation for
#' `a > 0`, so it is `NA` then. The paper's stochastic dynamic programme is
#' not implemented.
#'
#' @inheritParams rout2009
#' @param cost_ratio The cost ratio `R = C[s] / C[e]`, positive.
#' @param max_absent Largest number of absent surveys considered.
#'
#' @return A list with components `optimal_absent` (`d*`), `stop_time`
#'   (last sighting time plus `d*`), `approximation` (equation 5, or `NA`),
#'   `nec` (a data frame with columns `absent`, `time`, `posterior` and
#'   `nec`, the latter in units of `C[e]`), `n`, `sn`, `a`, `prior` and
#'   `cost_ratio`.
#'
#' @references
#' Rout, T. M., Salomon, Y., & McCarthy, M. A. (2009). Using sighting
#'   records to declare eradication of an invasive species. Journal of
#'   Applied Ecology, 46(1), 110-117.
#'
#' @export
rout2009_absences <- function(sd, cost_ratio, a = 0, prior = 0.5,
                              max_absent = 1000) {
  info <- rout2009_setup(sd, a, prior)
  if (missing(cost_ratio) || length(cost_ratio) != 1L || !is.numeric(cost_ratio) ||
      is.na(cost_ratio) || !is.finite(cost_ratio) || cost_ratio <= 0) {
    stop("`cost_ratio` must be a positive number.", call. = FALSE)
  }
  if (length(max_absent) != 1L || !is.numeric(max_absent) || is.na(max_absent) ||
      max_absent < 1 || max_absent != floor(max_absent)) {
    stop("`max_absent` must be a positive whole number.", call. = FALSE)
  }

  n <- info$n
  sn <- info$tn
  d <- 0:max_absent
  post <- vapply(d, function(k) {
    if (k == 0) return(1)
    b <- rout2009_bf(n, sn, sn + k, a)
    1 / (1 + (1 - prior) / (prior * b))
  }, numeric(1))
  nec <- (d - 1) * cost_ratio + post
  best <- d[which.min(nec)]

  approximation <- if (a == 0 && n >= 2) {
    sn * ((n - 1)^2 * prior / (sn * (1 - prior) * cost_ratio))^(1 / n) - sn
  } else {
    NA_real_
  }

  list(optimal_absent = best, stop_time = info$last_sight + best,
       approximation = approximation,
       nec = data.frame(absent = d, time = info$last_sight + d,
                        posterior = post, nec = nec),
       n = n, sn = sn, a = a, prior = prior, cost_ratio = cost_ratio)
}

#' Shared validation for rout2009() and rout2009_absences()
#' @keywords internal
#' @noRd
rout2009_setup <- function(sd, a, prior) {
  stopifnot(inherits(sd, "sighting_data"))
  if (length(prior) != 1L || !is.numeric(prior) || is.na(prior) ||
      prior <= 0 || prior >= 1) {
    stop("`prior` must be in (0, 1).", call. = FALSE)
  }
  if (length(a) != 1L || !is.numeric(a) || is.na(a) || a < 0 || a >= 1) {
    stop("`a` must be in [0, 1).", call. = FALSE)
  }
  if (!any(sd$count > 0)) stop("no sightings with count > 0.", call. = FALSE)
  origin <- min(sd$time)
  last_sight <- max(sd$time[sd$count > 0])
  if (last_sight <= origin) {
    stop("Rout et al. (2009) needs at least 1 sighting after the observation origin.",
         call. = FALSE)
  }
  list(origin = origin, last_sight = last_sight, tn = last_sight - origin,
       n = solow_event_count(sd))
}

#' Bayes factor of Rout et al. (2009) for times measured from the origin
#' @keywords internal
#' @noRd
rout2009_bf <- function(n, sn, s_total, a) {
  k <- n * (1 - a)
  ratio <- s_total / sn
  if (abs(k - 1) < 1e-10) 1 / log(ratio) else (k - 1) / (ratio^(k - 1) - 1)
}
