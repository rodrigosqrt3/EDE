#' Lee et al. (2014) Bayesian extinction probability with certain and
#' uncertain sightings
#'
#' Posterior probability that a species is extant at the end of `test_year`,
#' from a discrete record of certain sightings and, optionally, a second
#' record of uncertain sightings that may be false positives (Lee, McCarthy,
#' Wintle, Bode, Roberts & Burgman, 2014). The model is the OpenBUGS model of
#' the paper's Appendix S1, computed exactly instead of by MCMC.
#'
#' **Model** (Appendix S1; equations 5-12). The record is a series of
#' periods `i = 1, ..., T` (e.g. years). In each period at most one sighting
#' of each kind is recorded (equation 1). The extinction time is
#' `te ~ U(0, T)` given extinction (equation 12) and the species is present in
#' period `i` if it is extant at the end of the record or if `i <= te`
#' (`step(te - i)` in Appendix S1). A certain sighting occurs in a period with
#' probability `1 - exp(-mc)` while the species is present and `0` after
#' extinction (equation 5); an uncertain sighting occurs with probability
#' `1 - exp(-(mu + f))` while it is present and `1 - exp(-f)` after
#' extinction (equation 6). The two records are independent given `te`
#' (equations 7-8). The rates have independent uniform priors on
#' `(0, rate_max)` (equation 11; Appendix S1 uses `rate_max = 100` for `mc`,
#' `mu` and `f`), and the prior probability that the species is extant is
#' `prior` (`Extant ~ dbern(0.5)` in Appendix S1, `pi` in the paper). The
#' posterior is equation 2 with the Bayes factor of equation 3.
#'
#' **Exact computation.** Given extinction, the likelihood depends on `te`
#' only through the number of periods `k = floor(te)` with the species
#' present, which is uniform on `0, ..., T - 1`. For each `k`, the certain
#' likelihood is `(1/M) int_0^M (1 - e^-m)^S e^(-m (k - S)) dm` (equations
#' 9-10 in discrete form), an incomplete Beta function, and is zero if a
#' certain sighting occurs after period `k`. The uncertain likelihood is a
#' double integral over `mu` and `f`: the integral over `mu` is again an
#' incomplete Beta function and the integral over `f` is computed with
#' [stats::integrate()]. The likelihood under extinction (equation 10) is the
#' average over `k`; the likelihood under persistence is the case `k = T`.
#'
#' **Input.** Periods are the integer time units from the observation origin,
#' the earliest time in `certain` or `uncertain` (zero-count rows included),
#' to `test_year`, both included. A period is a sighting if its count is
#' positive (several sightings in one period count once). As in Appendix S1
#' (whose example record starts with a sighting), a sighting in the first
#' period counts as a sighting; add a zero-count row to start observation
#' before the first sighting. Uncertain sightings may occur anywhere in the
#' record, including before the first certain sighting.
#'
#' **Published conventions.** With these conventions the bars of Fig. 4 for
#' the paper's own model (Birdlife records, observation from `t0` in Table 3
#' to 2011) are reproduced. The bars of Fig. 3 and the 0.697 quoted for the
#' greater stick-nest rat are reproduced when the sighting at the origin is
#' not counted (the origin given as a zero-count row, as in Table 2) and the
#' record has `2011 - t0` periods, i.e. `test_year = 2010`.
#'
#' Package constructs: the exact computation (the paper uses MCMC), the
#' `extinction_time` output, and the requirement that `test_year` be later
#' than every sighting. The paper's eqn 14 (continuous extinction time) and
#' its prior for the false-positive rate built from the probability that a
#' sighting is true (eqn 13, Table 1; used for the light grey bars of Fig. 4)
#' are not implemented: the paper does not give enough detail to reproduce
#' those results. Appendix S1 has a single uncertain record; so has this
#' function (the uncertain sightings are grouped together, as in the dark
#' grey bars of Fig. 4).
#'
#' @param certain A [sighting_data] object with the certain sightings, with
#'   whole-number times. It may contain only zero-count rows (no certain
#'   sightings).
#' @param uncertain A [sighting_data] object with the uncertain sightings, or
#'   `NULL` (certain sightings only).
#' @param test_year Last period of observation. Must be a whole number later
#'   than every sighting.
#' @param prior Prior probability that the species is extant at the end of
#'   `test_year`, in (0, 1).
#' @param rate_max Upper bound `M` of the uniform priors on the sighting
#'   rates (Appendix S1: 100).
#'
#' @return A list of class `ede_bayes` with components `bayes_factor`
#'   (`l(s | H0) / l(s | H1)`, equation 3), `posterior` (probability that the
#'   species is extant, equation 2), `prior`, `test_year`, `method`,
#'   `periods` (`T`), `rate_max`, and `extinction_time`: a data frame with,
#'   for each period after the last certain sighting, the posterior
#'   probability that it is the first period in which the species is absent
#'   (these probabilities sum to `1 - posterior`).
#'
#' @references
#' Lee, T. E., McCarthy, M. A., Wintle, B. A., Bode, M., Roberts, D. L. &
#'   Burgman, M. A. (2014). Inferring extinctions from sighting records of
#'   variable reliability. Journal of Applied Ecology, 51(1), 251-258.
#'
#' @examples
#' # Greater stick-nest rat (Lee et al. 2014, Table 2), observation from 1835
#' yrs <- c(1844, 1857, 1879, 1896, 1897, 1907, 1920, 1921, 1930, 1933,
#'          1938, 1986)
#' rat <- sighting_data(data.frame(year = c(1835, yrs),
#'                                 count = c(0, rep(1, 12))))
#' lee2014(rat, test_year = 2010)$posterior  # Fig. 3: 0.697
#' # the 1986 sighting treated as uncertain
#' rat_c <- sighting_data(data.frame(year = c(1835, yrs[-12]),
#'                                   count = c(0, rep(1, 11))))
#' rat_u <- sighting_data(data.frame(year = 1986, count = 1))
#' lee2014(rat_c, rat_u, test_year = 2010)
#'
#' @export
lee2014 <- function(certain, uncertain = NULL, test_year, prior = 0.5,
                    rate_max = 100) {
  if (!inherits(certain, "sighting_data")) {
    stop("`certain` must be a sighting_data object.", call. = FALSE)
  }
  if (!is.null(uncertain) && !inherits(uncertain, "sighting_data")) {
    stop("`uncertain` must be a sighting_data object or NULL.", call. = FALSE)
  }
  if (length(prior) != 1L || !is.numeric(prior) || is.na(prior) ||
      prior <= 0 || prior >= 1) {
    stop("`prior` must be in (0, 1).", call. = FALSE)
  }
  if (length(rate_max) != 1L || !is.numeric(rate_max) || !is.finite(rate_max) ||
      rate_max <= 0) {
    stop("`rate_max` must be a positive number.", call. = FALSE)
  }
  times <- c(certain$time, if (!is.null(uncertain)) uncertain$time)
  if (any(times != round(times))) {
    stop("lee2014() requires whole-number times (one row per period).",
         call. = FALSE)
  }
  pos <- c(certain$time[certain$count > 0],
           if (!is.null(uncertain)) uncertain$time[uncertain$count > 0])
  if (length(pos) == 0L) stop("no sightings with count > 0.", call. = FALSE)
  check_test_year(test_year, max(pos))
  if (test_year != round(test_year)) {
    stop("`test_year` must be a whole number.", call. = FALSE)
  }
  origin <- min(times)
  T <- test_year - origin + 1
  yc <- seq_len(T) %in% (certain$time[certain$count > 0] - origin + 1)
  yu <- if (is.null(uncertain)) NULL else
    seq_len(T) %in% (uncertain$time[uncertain$count > 0] - origin + 1)

  S <- sum(yc)
  last <- if (S > 0) max(which(yc)) else 0
  cum_u <- if (is.null(yu)) NULL else c(0, cumsum(yu))
  nu <- sum(yu)
  loglik <- function(k) {
    v <- lee2014_lcert(S, k, rate_max)
    if (!is.null(yu)) {
      a <- cum_u[k + 1]
      cc <- nu - a
      v <- v + lee2014_lunc(a, k - a, cc, T - k - cc, rate_max)
    }
    v
  }
  ks <- seq(last, T - 1)  # last < T: test_year is later than every sighting
  l0 <- loglik(T)
  l1 <- vapply(ks, loglik, numeric(1))
  mx <- max(c(l0, l1))
  w1 <- exp(l1 - mx) / T          # l(s, k | H1) for each k, scaled
  bf <- exp(l0 - mx) / sum(w1)
  posterior <- prior * bf / (prior * bf + 1 - prior)
  probability <- (1 - posterior) * w1 / sum(w1)
  method <- if (is.null(uncertain)) {
    "Lee et al. (2014) Bayesian model, certain sightings"
  } else {
    "Lee et al. (2014) Bayesian model, certain and uncertain sightings"
  }
  structure(
    list(bayes_factor = bf, posterior = posterior, prior = prior,
         test_year = test_year, method = method, periods = T,
         rate_max = rate_max,
         extinction_time = data.frame(time = origin + ks,
                                      probability = probability)),
    class = "ede_bayes"
  )
}

#' log of int_{x0}^{x1} (1 - x)^a x^(b - 1) dx, from log(x0) and log(x1)
#'
#' Vectorised over `lx0` and `lx1` (`a`, `b` scalar). With `x = exp(-u)` this
#' is int (1 - exp(-u))^a exp(-u b) du over `u` in (-lx1, -lx0).
#' @keywords internal
#' @noRd
lee2014_lbint <- function(a, b, lx0, lx1) {
  if (b > 0) {
    l1 <- stats::pbeta(exp(lx1), b, a + 1, log.p = TRUE)
    l0 <- stats::pbeta(exp(lx0), b, a + 1, log.p = TRUE)
    return(lbeta(b, a + 1) + l1 + log1p(-exp(l0 - l1)))
  }
  # b = 0: int (1 - x)^a / x dx
  #      = log(x1 / x0) - sum_{j=1}^a ((1 - x0)^j - (1 - x1)^j) / j
  j <- seq_len(a)
  val <- vapply(seq_along(lx0), function(i) {
    s <- if (a > 0) sum(((1 - exp(lx0[i]))^j - (1 - exp(lx1[i]))^j) / j) else 0
    lx1[i] - lx0[i] - s
  }, numeric(1))
  bad <- !(val > 1e-6 * (lx1 - lx0))
  for (i in which(bad)) {
    # cancellation: integrate (1 - exp(-u))^a over u directly
    g <- function(u) exp(a * log1p(-exp(-u)))
    val[i] <- stats::integrate(g, -lx1[i], -lx0[i], rel.tol = 1e-10)$value
  }
  log(val)
}

#' log of (1/M) int_0^M (1 - e^-m)^S e^(-m (k - S)) dm; -Inf if k < S
#' @keywords internal
#' @noRd
lee2014_lcert <- function(S, k, M) {
  if (k < S) return(-Inf)
  lee2014_lbint(S, k - S, -M, 0) - log(M)
}

#' log of (1/M^2) int_0^M int_0^M (1 - e^-(m+f))^a e^(-(m+f) b)
#'   (1 - e^-f)^c e^(-f d) dm df
#'
#' `a`, `b`: uncertain sightings and non-sightings while present; `c`, `d`:
#' after extinction. The integral over `m` is closed form (lee2014_lbint).
#' @keywords internal
#' @noRd
lee2014_lunc <- function(a, b, c, d, M) {
  logh <- function(f) {
    lc <- if (c > 0) c * log1p(-exp(-f)) else 0
    lc - f * d + lee2014_lbint(a, b, -f - M, -f)
  }
  br <- c(0, M * 10^seq(-8, 0, by = 0.25))
  mids <- sqrt(br[-c(1, length(br))] * br[-(1:2)])
  mx <- max(logh(c(br[-1], mids)))
  tot <- 0
  for (i in seq_len(length(br) - 1L)) {
    tot <- tot + stats::integrate(function(f) exp(logh(f) - mx), br[i],
                                  br[i + 1L], rel.tol = 1e-10)$value
  }
  mx + log(tot) - 2 * log(M)
}
