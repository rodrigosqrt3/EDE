#' Solow et al. (2012) Bayes factor with certain and uncertain sightings
#'
#' Not exported: the published Bayes factors of this paper and of Solow &
#' Beet (2014) are not yet reproduced (see "Reproduction status" below).
#'
#' Bayesian test of persistence for a record that mixes certain sightings
#' (known to be valid) with uncertain sightings (valid or invalid), following
#' Solow, Smith, Burgman, Rout, Wintle & Roberts (2012). The record is assumed
#' to start with a period `(0, tL)` that contains only certain sightings, `tL`
#' being the last certain sighting. Valid sightings (certain, and uncertain
#' ones before extinction) follow a stationary Poisson process with rate
#' `lambda`; invalid sightings follow an independent Poisson process with rate
#' `lambda * (1 - omega) / omega` over `(tL, T)`, so `omega` is the expected
#' proportion of valid uncertain sightings before extinction. Conditional on
#' extinction at `tau` (`tL <= tau <= T`), the likelihood is equation 5 of the
#' paper,
#' `lambda^n omega^(m - n) (1 - omega)^m(tau) exp(-lambda (tau + (1 - omega) / omega (T - tL)))`,
#' with `n` the total number of sightings, `m` the number of certain sightings
#' and `m(tau)` the number of sightings after `tau`. With the noninformative
#' prior `p(lambda, omega) = 1 / lambda` (equation 6), `lambda` is integrated
#' out analytically and `omega` numerically (equation 7). The likelihood
#' under extinction (equation 3) integrates equation 7 against the prior of
#' the extinction time over `(tL, T)`; the likelihood under persistence is
#' equation 7 evaluated at `tau = T`. Their ratio is the Bayes factor in
#' favour of extinction (equation 2).
#'
#' Times are measured from the observation origin, the earliest time in
#' `certain`. As in [solow1993()], one certain sighting at the origin defines
#' it and is omitted; every other sighting event is counted (a time unit with
#' several sightings contributes several events). The test year is the end
#' `T` of the observation period. Uncertain sightings must all be later than
#' the last certain sighting (the paper assumes uncertain sightings before it
#' were omitted from the record); use [solow_beet2014_bayes()] when the two
#' kinds of sightings overlap in time.
#'
#' Package constructs: the paper reports the Bayes factor in favour of
#' extinction (returned as `bayes_factor_extinction`); for consistency with
#' [solow1993_bayes()], `bayes_factor` is its reciprocal, the Bayes factor in
#' favour of persistence, and `posterior` is the posterior probability of
#' persistence for prior probability `prior`. The integrals over `tau` are
#' computed with [stats::integrate()] between consecutive uncertain sighting
#' times (where the likelihood is discontinuous), and the integral over
#' `omega` on the log scale.
#'
#' Reproduction status: the posterior of `omega` reproduces the paper's mode
#' (0.67) and quartiles (0.64, 0.77 versus 0.64, 0.78 here) for the
#' Ivory-billed Woodpecker, but the published Bayes factors (4.03 for the
#' truncated exponential prior, 4.97 for the uniform prior) are not
#' reproduced exactly (3.28 and 3.87 here); the paper does not state the
#' time conventions behind its numerical evaluation, to which the Bayes
#' factor is sensitive.
#'
#' @param certain A [sighting_data] object with the certain sightings. Its
#'   earliest time defines the observation origin.
#' @param uncertain A [sighting_data] object with the uncertain sightings, or
#'   `NULL` if there are none. Zero-count rows are ignored.
#' @param test_year End `T` of the observation period. Must be later than
#'   every sighting.
#' @param prior Prior probability that the species persisted to `test_year`,
#'   in (0, 1).
#' @param tau_prior Prior of the extinction time given extinction during the
#'   observation period: `"uniform"` on `(origin, test_year)`,
#'   `"exponential"` (the truncated exponential of equation 8, set through
#'   `exp_quantile`), or a function of time returning a non-negative density,
#'   normalised numerically over `(origin, test_year)`.
#' @param exp_quantile For `tau_prior = "exponential"`: a numeric vector
#'   `c(time, prob)` giving the prior probability `prob` that extinction
#'   occurred before `time`, conditional on extinction during the observation
#'   period. Solow et al. (2012) used `c(1940, 0.5)`. A `prob` below
#'   `(time - origin) / (test_year - origin)` gives an increasing density.
#' @param data_out `FALSE` (default) for the Bayes factor; `TRUE` or
#'   `"extinction"` for the posterior density of the extinction time given
#'   extinction (equation 9); `"omega"` for the posterior density of `omega`
#'   given extinction (equations 10-11).
#'
#' @return A list of class `ede_bayes` with components `bayes_factor` (in
#'   favour of persistence), `posterior` (of persistence), `prior`,
#'   `test_year`, `method`, plus `bayes_factor_extinction` (the paper's
#'   `B(t)`) and `tau_prior`. With `data_out = "extinction"`, a data frame
#'   with columns `time`, `likelihood_ratio` (`p(t | tau) / p(t | persistence)`)
#'   and `density`; with `data_out = "omega"`, a data frame with columns
#'   `omega` and `density`.
#'
#' @references
#' Solow, A., Smith, W., Burgman, M., Rout, T., Wintle, B. & Roberts, D.
#'   (2012). Uncertain sightings and the extinction of the Ivory-billed
#'   Woodpecker. Conservation Biology, 26(1), 180-184.
#'
#' @seealso [solow_beet2014_bayes()] for records in which certain and
#'   uncertain sightings overlap in time.
#' @keywords internal
#' @noRd
solow2012_bayes <- function(certain, uncertain = NULL, test_year, prior = 0.5,
                            tau_prior = c("uniform", "exponential"),
                            exp_quantile = NULL, data_out = FALSE) {
  rec <- uncertain_record(certain, uncertain, test_year, "Solow et al. (2012)")
  if (rec$nu > 0 && min(rec$u_time) <= rec$tL) {
    stop("Solow et al. (2012) requires every uncertain sighting to be later ",
         "than the last certain sighting; use solow_beet2014_bayes() for ",
         "overlapping records.", call. = FALSE)
  }
  model <- list(a = rec$nc, n = rec$nc + rec$nu, L = rec$T - rec$tL,
                certain_power = 0)
  uncertain_bayes(rec, model, prior, tau_prior, exp_quantile, data_out,
                  "Solow et al. (2012) Bayes factor")
}

#' Solow & Beet (2014) Bayes factors with overlapping uncertain sightings
#'
#' Extension of [solow2012_bayes()] by Solow & Beet (2014) in which certain
#' and uncertain sightings may occur anywhere in the observation period
#' `(0, T)`. Two models are given. In both, conditional on extinction at
#' `tau`, valid uncertain sightings occur before `tau` and invalid sightings
#' follow a stationary Poisson process with rate `lambda (1 - omega) / omega`
#' over the whole period `(0, T)`.
#'
#' * Model 1: valid sightings, certain or uncertain, come from the same
#'   Poisson process with rate `lambda` over `(0, tau)`. The likelihood is
#'   equation 6,
#'   `lambda^n exp(-lambda (tau + (1 - omega) / omega T)) omega^(-nu) (1 - omega)^(n - n(tau))`,
#'   with `nu` the number of uncertain sightings and `n - n(tau)` the number
#'   of sightings after `tau`. Integrating `lambda` against `1 / lambda`
#'   gives equation 7, and `omega` is integrated numerically (equation 8).
#' * Model 2: certain sightings follow an independent Poisson process with
#'   rate `M` (prior `1 / M`) over `(0, tau)`, so
#'   `p(t | tau) = p(tc | tau) p(tu | tau)` (equation 9) with
#'   `p(tc | tau) = (nc - 1)! / tau^nc` (equation 10) and `p(tu | tau)` given
#'   by equation 11, the model-1 integral applied to the uncertain sightings
#'   alone.
#'
#' The Bayes factor in favour of extinction is formed exactly as in Solow
#' et al. (2012): equations 2 and 3, integrating over `tau` from the last
#' certain sighting `tL` to `T`. With no uncertain sightings, model 2
#' reduces to the Bayes factor of [solow1993_bayes()] under a uniform prior.
#' Times, the origin, event counts, `prior`, `tau_prior` and the returned
#' quantities are as in [solow2012_bayes()].
#'
#' Reproduction status: for the Ivory-billed Woodpecker the paper reports
#' Bayes factors in favour of extinction of 0.13, 0.19 and 0.007 (model 1)
#' and about 4e6, 4e6 and 3e6 (model 2) for the truncated exponential prior
#' with probability 0.5 before 1940, the uniform prior and the truncated
#' exponential prior with probability 0.9 before 1940. This implementation
#' gives 0.24, 0.36 and 0.013 (model 1) and 3.3e6, 3.1e6 and 1.8e6
#' (model 2): the same conclusions, but not the same values; the paper does
#' not state the time conventions behind its numerical evaluation.
#'
#' @inheritParams solow2012_bayes
#' @param model `1` (certain and valid uncertain sightings from the same
#'   process) or `2` (independent processes).
#'
#' @return As [solow2012_bayes()].
#'
#' @references
#' Solow, A. R. & Beet, A. R. (2014). On uncertain sightings and inference
#'   about extinction. Conservation Biology, 28(4). doi:10.1111/cobi.12309
#'
#' Solow, A., Smith, W., Burgman, M., Rout, T., Wintle, B. & Roberts, D.
#'   (2012). Uncertain sightings and the extinction of the Ivory-billed
#'   Woodpecker. Conservation Biology, 26(1), 180-184.
#'
#' @keywords internal
#' @noRd
solow_beet2014_bayes <- function(certain, uncertain = NULL, test_year, model,
                                 prior = 0.5,
                                 tau_prior = c("uniform", "exponential"),
                                 exp_quantile = NULL, data_out = FALSE) {
  if (missing(model) || length(model) != 1L || !(model %in% c(1, 2))) {
    stop("`model` must be 1 or 2.", call. = FALSE)
  }
  rec <- uncertain_record(certain, uncertain, test_year, "Solow & Beet (2014)")
  spec <- if (model == 1) {
    list(a = rec$nc, n = rec$nc + rec$nu, L = rec$T, certain_power = 0)
  } else {
    list(a = 0, n = rec$nu, L = rec$T, certain_power = rec$nc)
  }
  uncertain_bayes(rec, spec, prior, tau_prior, exp_quantile, data_out,
                  sprintf("Solow & Beet (2014) model %d Bayes factor", model))
}

# ---------------------------------------------------------------------------
# Internal engine
#
# All three models share the conditional log-likelihood (up to constants)
#   log p(t | tau) = -cp * log(tau) + log int_0^1 exp(f(omega, tau)) d omega,
#   f = a log(omega) + k(tau) log(1 - omega) - n log(L + (tau - L) omega),
# obtained from tau + (1 - omega) / omega * L = (L + (tau - L) omega) / omega:
#   Solow et al. (2012) eq 7:  a = m, n = m + nu, L = T - tL, cp = 0
#   Solow & Beet model 1 eq 7: a = nc, n = nc + nu, L = T, cp = 0
#   Solow & Beet model 2 eqs 10-11: a = 0, n = nu, L = T, cp = nc
# with k(tau) the number of uncertain sightings after tau.
# ---------------------------------------------------------------------------

#' Validate a certain/uncertain record and convert it to relative times
#' @keywords internal
#' @noRd
uncertain_record <- function(certain, uncertain, test_year, label) {
  if (!inherits(certain, "sighting_data")) {
    stop("`certain` must be a sighting_data object.", call. = FALSE)
  }
  if (!is.null(uncertain) && !inherits(uncertain, "sighting_data")) {
    stop("`uncertain` must be a sighting_data object or NULL.", call. = FALSE)
  }
  if (!any(certain$count > 0)) stop("no certain sightings with count > 0.", call. = FALSE)
  origin <- min(certain$time)
  last_certain <- max(certain$time[certain$count > 0])
  nc <- solow_event_count(certain)
  if (last_certain <= origin || nc < 1) {
    stop(label, " needs at least 1 certain sighting after the observation origin.",
         call. = FALSE)
  }
  u_time <- numeric(0)
  u_count <- numeric(0)
  if (!is.null(uncertain)) {
    keep <- uncertain$count > 0
    u_time <- uncertain$time[keep]
    u_count <- uncertain$count[keep]
    if (any(u_time <= origin)) {
      stop("uncertain sightings must be later than the observation origin ",
           "(the earliest time in `certain`).", call. = FALSE)
    }
  }
  check_test_year(test_year, max(c(last_certain, u_time)))
  list(origin = origin, nc = nc, nu = sum(u_count), tL = last_certain - origin,
       T = test_year - origin, u_time = u_time - origin, u_count = u_count,
       test_year = test_year)
}

#' Prior density of the extinction time on (0, T), relative time
#' @keywords internal
#' @noRd
uncertain_tau_prior <- function(tau_prior, exp_quantile, rec) {
  Tt <- rec$T
  if (is.function(tau_prior)) {
    f <- function(x) {
      v <- tau_prior(x + rec$origin)
      if (length(v) != length(x) || any(!is.finite(v)) || any(v < 0)) {
        stop("`tau_prior` must return finite non-negative densities.", call. = FALSE)
      }
      v
    }
    z <- stats::integrate(f, 0, Tt, rel.tol = 1e-10, subdivisions = 1000L)$value
    if (!(z > 0)) stop("`tau_prior` integrates to zero.", call. = FALSE)
    return(list(density = function(x) f(x) / z, label = "user-supplied"))
  }
  tau_prior <- match.arg(tau_prior, c("uniform", "exponential"))
  if (tau_prior == "uniform") {
    return(list(density = function(x) rep(1 / Tt, length(x)), label = "uniform"))
  }
  if (is.null(exp_quantile) || !is.numeric(exp_quantile) ||
      length(exp_quantile) != 2L || anyNA(exp_quantile)) {
    stop("`exp_quantile` must be c(time, prob) for the exponential prior.",
         call. = FALSE)
  }
  a <- exp_quantile[1] - rec$origin
  p <- exp_quantile[2]
  if (a <= 0 || a >= Tt) {
    stop("`exp_quantile` time must lie inside the observation period.", call. = FALSE)
  }
  if (p <= 0 || p >= 1) stop("`exp_quantile` prob must be in (0, 1).", call. = FALSE)
  label <- sprintf("truncated exponential, P(extinct before %s) = %g",
                   format(exp_quantile[1]), p)
  # P(tau < a) = (1 - exp(-g a)) / (1 - exp(-g T)), increasing in g (eq 8)
  cdf_a <- function(g) {
    if (abs(g * Tt) < 1e-8) return(a / Tt)
    expm1(-g * a) / expm1(-g * Tt)
  }
  if (abs(p - a / Tt) < 1e-12) {
    return(list(density = function(x) rep(1 / Tt, length(x)), label = label))
  }
  bound <- 600 / Tt
  if ((cdf_a(-bound) - p) * (cdf_a(bound) - p) > 0) {
    stop("`exp_quantile` is too extreme for a truncated exponential prior.",
         call. = FALSE)
  }
  g <- stats::uniroot(function(g) cdf_a(g) - p, c(-bound, bound),
                      tol = 1e-14)$root
  density <- function(x) {
    # g exp(-g x) / (1 - exp(-g T)), written to avoid overflow for g < 0
    if (g > 0) g * exp(-g * x) / (-expm1(-g * Tt))
    else g * exp(-g * (x - Tt)) / expm1(g * Tt)
  }
  list(density = density, label = label)
}

#' a * log(x) with the convention 0 * log(0) = 0
#' @keywords internal
#' @noRd
uncertain_xlog <- function(a, x) if (a == 0) 0 * x else a * log(x)

#' log of the integrand over omega (see header), excluding -cp log(tau)
#' @keywords internal
#' @noRd
uncertain_logf <- function(w, tau, k, spec) {
  uncertain_xlog(spec$a, w) + uncertain_xlog(k, 1 - w) -
    uncertain_xlog(spec$n, spec$L + (tau - spec$L) * w)
}

#' log p(t | tau) up to a constant, for a single tau
#' @keywords internal
#' @noRd
uncertain_loglik <- function(tau, rec, spec) {
  k <- sum(rec$u_count[rec$u_time > tau])
  certain_part <- -spec$certain_power * log(tau)
  if (spec$a == 0 && k == 0 && spec$n == 0) return(certain_part)
  f <- function(w) uncertain_logf(w, tau, k, spec)
  # locate the global maximum on a grid, then refine it
  grid <- seq(0.0025, 0.9975, by = 0.0025)
  fg <- f(grid)
  i <- which.max(fg)
  lower <- if (i == 1L) 0 else grid[i - 1L]
  upper <- if (i == length(grid)) 1 else grid[i + 1L]
  opt <- stats::optimize(f, c(lower, upper), maximum = TRUE, tol = 1e-10)
  fmax <- max(opt$objective, fg[i])
  wmax <- if (opt$objective >= fg[i]) opt$maximum else grid[i]
  g <- function(w) exp(f(w) - fmax)
  v <- 0
  if (wmax > 0) {
    v <- v + stats::integrate(g, 0, wmax, rel.tol = 1e-10,
                              subdivisions = 1000L)$value
  }
  if (wmax < 1) {
    v <- v + stats::integrate(g, wmax, 1, rel.tol = 1e-10,
                              subdivisions = 1000L)$value
  }
  certain_part + fmax + log(v)
}

#' Shared driver: Bayes factor, posterior and data_out curves
#' @keywords internal
#' @noRd
uncertain_bayes <- function(rec, spec, prior, tau_prior, exp_quantile,
                            data_out, method) {
  if (length(prior) != 1L || !is.numeric(prior) || is.na(prior) ||
      prior <= 0 || prior >= 1) {
    stop("`prior` must be in (0, 1).", call. = FALSE)
  }
  if (isTRUE(data_out)) data_out <- "extinction"
  if (!identical(data_out, FALSE) &&
      !(is.character(data_out) && length(data_out) == 1L &&
        data_out %in% c("extinction", "omega"))) {
    stop("`data_out` must be FALSE, TRUE, \"extinction\" or \"omega\".",
         call. = FALSE)
  }
  pri <- uncertain_tau_prior(tau_prior, exp_quantile, rec)

  l0 <- uncertain_loglik(rec$T, rec, spec)
  ratio <- function(x) {
    vapply(x, function(s) exp(uncertain_loglik(s, rec, spec) - l0), numeric(1))
  }
  # the likelihood is a step function of k(tau): integrate between jumps
  breaks <- sort(unique(c(rec$tL, rec$u_time[rec$u_time > rec$tL], rec$T)))
  bf_ext <- 0
  for (i in seq_len(length(breaks) - 1L)) {
    bf_ext <- bf_ext + stats::integrate(
      function(x) ratio(x) * pri$density(x), breaks[i], breaks[i + 1L],
      rel.tol = 1e-8, subdivisions = 1000L
    )$value
  }

  if (identical(data_out, "extinction")) {
    tau <- seq(rec$tL, rec$T, length.out = 501L)
    lr <- ratio(tau)
    return(data.frame(time = rec$origin + tau, likelihood_ratio = lr,
                      density = lr * pri$density(tau) / bf_ext))
  }
  if (identical(data_out, "omega")) {
    omega <- seq(0.005, 0.995, by = 0.005)
    dens <- vapply(omega, function(w) {
      tot <- 0
      for (i in seq_len(length(breaks) - 1L)) {
        k <- sum(rec$u_count[rec$u_time > breaks[i]])
        h <- function(x) {
          exp(uncertain_logf(w, x, k, spec) - spec$certain_power * log(x) - l0) *
            pri$density(x)
        }
        tot <- tot + stats::integrate(h, breaks[i], breaks[i + 1L],
                                      rel.tol = 1e-8)$value
      }
      tot
    }, numeric(1))
    return(data.frame(omega = omega, density = dens / bf_ext))
  }

  b <- 1 / bf_ext
  structure(
    list(bayes_factor = b,
         posterior = 1 / (1 + (1 - prior) / (prior * b)),
         prior = prior, test_year = rec$test_year, method = method,
         bayes_factor_extinction = bf_ext, tau_prior = pri$label),
    class = "ede_bayes"
  )
}
