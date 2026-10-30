#' Kodikara et al. (2020) Bayesian extinction year from yearly sighting records
#'
#' Posterior distribution of the extinction year from a record of yearly
#' certain and (optionally) uncertain sightings, following Kodikara,
#' Demirhan, Wang, Solow & Stone (2020). Each year `t = 1, ..., T` of the
#' observation period is a trial. `tau_E` is the first year following
#' extinction: the species is extant in years `1, ..., tau_E - 1`, so
#' `tau_E` is at least one year after the last certain sighting `c_n`.
#'
#' * **Model 1** (certain sightings only, section 2.1). An extant species is
#'   sighted (with certainty) in a year with probability `pc`. With `nc`
#'   certain sighting years, the likelihood is equation 3:
#'   `0` for `tau_E <= c_n`, `pc^nc (1 - pc)^(tau_E - 1 - nc)` for
#'   `c_n < tau_E <= T`, and `pc^nc (1 - pc)^(T - nc)` for `tau_E > T`.
#' * **Model 2** (certain and uncertain sightings, section 2.2). Before
#'   extinction, valid and invalid uncertain sightings occur independently
#'   with probabilities `puv` and `pui`, so an uncertain sighting occurs with
#'   probability `pu = puv (1 - pui) + pui (1 - puv) + puv pui`, and is
#'   recorded only in years without a certain sighting (the "censoring"
#'   process). After extinction, only invalid uncertain sightings occur, with
#'   probability `pui`. With `nu(tau_E)` uncertain sightings before `tau_E`
#'   and `Nu` in total, the likelihood is equation 9:
#'   `pc^nc ((1 - pc) pu)^nu(tau_E) ((1 - pc)(1 - pu))^(tau_E - 1 - nc - nu(tau_E))`
#'   `pui^(Nu - nu(tau_E)) (1 - pui)^(T - (tau_E - 1) - (Nu - nu(tau_E)))` for
#'   `c_n < tau_E <= T`, and
#'   `pc^nc ((1 - pc) pu)^Nu ((1 - pc)(1 - pu))^(T - nc - Nu)` for
#'   `tau_E > T`.
#' * **Model 3** (section 4.2, Table 4): Model 2 with `puv = pc`.
#'
#' The prior of `tau_E` is geometric with parameter `theta` (equation 4),
#' `p(tau_E | theta) = (1 - theta)^(tau_E - 1) theta`, and `theta`, `pc`,
#' `puv` and `pui` are standard uniform (equations 5-6); the joint posterior
#' is equation 7.
#'
#' The paper samples this posterior with JAGS. Because `tau_E` is discrete
#' and all priors are uniform, the posterior of `tau_E` is computed here
#' exactly instead (a package construct): integrating `theta` gives the
#' marginal prior `p(tau_E) = 1 / (tau_E (tau_E + 1))`, whose tail beyond
#' `T` has mass `1 / (T + 1)`; the integrals over `pc`, `puv` and `pui` are
#' finite sums of Beta functions, obtained by writing
#' `1 - (1 - puv)(1 - pui) = puv + (1 - puv) pui` and expanding its power
#' binomially (all terms positive). Since the likelihood is constant for
#' `tau_E > T`, the posterior of `tau_E` beyond `T` is proportional to
#' `1 / (tau_E (tau_E + 1))`, so its quantiles there are analytic. The
#' posterior median and highest posterior density interval (HDI, the
#' shortest run of consecutive years with posterior probability at least
#' `level`) of `tau_E` are reported as years, as in Tables 2-4 of the paper.
#'
#' **Input.** Periods are the integer years from the earliest time in
#' `certain` to `test_year` (`T`); include a zero-count row in `certain` when
#' observation began before the first certain sighting (as for the O'ahu
#' 'Alauahio record, which starts in 1837). A year is a sighting year if its
#' count is positive (several sightings in a year count once); unlisted years
#' have no sighting. Uncertain sightings in a year that also has a certain
#' sighting are ignored, as the censoring process of Model 2 assumes they
#' cannot be recorded.
#'
#' @param certain A [sighting_data] object with whole-number times (years)
#'   and the certain sightings. Its earliest time is the first year of
#'   observation. At least one certain sighting is required.
#' @param uncertain A [sighting_data] object with the uncertain sightings,
#'   or `NULL`. Ignored by Model 1.
#' @param test_year Last year `T` of the observation period. Must be later
#'   than the last certain sighting and not earlier than the last uncertain
#'   sighting.
#' @param model `1`, `2` or `3` (see Description). Defaults to `1` without
#'   uncertain sightings and `2` otherwise.
#' @param level Probability content of the HDIs.
#'
#' @return A list of class `ede_bayes` with components:
#'
#' * `posterior`: posterior probability `p(tau_E > T | S)` that the species
#'   was extant in every year of the record (one minus the paper's
#'   `p(tau_E <= T | S)`).
#' * `prior`: the implied prior of that event, `1 / (T + 1)`.
#' * `bayes_factor`: posterior odds over prior odds of persistence.
#' * `test_year`, `method` and `model`.
#' * `extinction_year`: posterior median and `level` HDI of `tau_E`, as
#'   years (`c(median, lower, upper)`); `level`.
#' * `extinction_time`: data frame with the posterior probability of each
#'   value of `tau_E` (`time`, as a year) from the year after the last
#'   certain sighting to `test_year`.
#' * `p_certain`: posterior median and `level` HDI of `pc` (Models 1 and 2;
#'   `NA` for Model 3).
#'
#' @references
#' Kodikara, S., Demirhan, H., Wang, Y., Solow, A. & Stone, L. (2020).
#'   Inferring extinction year using a Bayesian approach. Methods in Ecology
#'   and Evolution, 11(8), 964-973. doi:10.1111/2041-210X.13408
#'
#' @examples
#' # Ivory-billed Woodpecker, physical-evidence sightings (Kodikara et al.,
#' # 2020, Table 2)
#' ibw_c <- sighting_data(data.frame(
#'   year = c(1898:1902, 1904:1910, 1913, 1914, 1917, 1924, 1925, 1932,
#'            1935, 1938, 1939),
#'   count = 1))
#' fit <- kodikara2020(ibw_c, test_year = 2010)
#' fit$extinction_year
#'
#' @export
kodikara2020 <- function(certain, uncertain = NULL, test_year, model = NULL,
                         level = 0.95) {
  if (!inherits(certain, "sighting_data")) {
    stop("`certain` must be a sighting_data object.", call. = FALSE)
  }
  if (!is.null(uncertain) && !inherits(uncertain, "sighting_data")) {
    stop("`uncertain` must be a sighting_data object or NULL.", call. = FALSE)
  }
  if (is.null(model)) model <- if (is.null(uncertain)) 1 else 2
  if (length(model) != 1L || !is.numeric(model) || !(model %in% 1:3)) {
    stop("`model` must be 1, 2 or 3.", call. = FALSE)
  }
  if (length(level) != 1L || !is.numeric(level) || is.na(level) ||
      level <= 0 || level >= 1) {
    stop("`level` must be in (0, 1).", call. = FALSE)
  }
  if (!any(certain$count > 0)) {
    stop("no certain sightings with count > 0.", call. = FALSE)
  }
  whole <- function(x) all(x == round(x))
  if (!whole(certain$time) || (!is.null(uncertain) && !whole(uncertain$time))) {
    stop("kodikara2020() requires whole-number times (years).", call. = FALSE)
  }
  origin <- min(certain$time)
  last_certain <- max(certain$time[certain$count > 0])
  check_test_year(test_year, last_certain)
  if (test_year != round(test_year)) {
    stop("`test_year` must be a whole number.", call. = FALSE)
  }
  yc <- as.integer(expand_record(certain, test_year)$count > 0)
  T <- length(yc)
  yu <- integer(T)
  if (model > 1) {
    if (is.null(uncertain)) {
      stop("Models 2 and 3 need `uncertain` sightings.", call. = FALSE)
    }
    ut <- uncertain$time[uncertain$count > 0]
    if (any(ut < origin) || any(ut > test_year)) {
      stop("uncertain sightings must lie between the earliest time in ",
           "`certain` and `test_year`.", call. = FALSE)
    }
    yu[ut - origin + 1] <- 1L
    yu[yc == 1L] <- 0L                          # censoring: certain years
  }

  fit <- kodikara2020_posterior(yc, yu, model)
  ext <- kodikara2020_summary(fit$w, fit$k, T, level)
  pers <- fit$w[length(fit$w)]
  prior <- 1 / (T + 1)
  bf <- if (pers >= 1) Inf else (pers / (1 - pers)) / (prior / (1 - prior))
  finite <- seq_len(length(fit$k) - 1L)
  labels <- c("certain sightings only", "certain and uncertain sightings",
              "certain and uncertain sightings, puv = pc")
  structure(
    list(bayes_factor = bf, posterior = pers, prior = prior,
         test_year = test_year,
         method = sprintf("Kodikara et al. (2020) Model %d (%s)", model,
                          labels[model]),
         model = model,
         extinction_year = origin - 1 + ext, level = level,
         extinction_time = data.frame(time = origin - 1 + fit$k[finite],
                                      probability = fit$w[finite]),
         p_certain = if (model < 3) {
           kodikara2020_pc(fit$w, fit$pc_a, fit$pc_b, level)
         } else {
           c(median = NA_real_, lower = NA_real_, upper = NA_real_)
         }),
    class = "ede_bayes"
  )
}

#' log of sum(exp(x))
#' @keywords internal
#' @noRd
kodikara2020_lse <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

#' Posterior weights of tau_E = tn + 1, ..., T and of tau_E > T (last entry)
#'
#' For index k1 = min(tau_E, T + 1) - 1 years before extinction, with
#' a = nu(tau_E) uncertain sightings and b = k1 - nc - a empty years before
#' extinction, c = Nu - a uncertain and d = T - k1 - c empty years after:
#' Model 1: B(nc + 1, k1 - nc + 1).
#' Model 2: B(nc + 1, a + b + 1) *
#'   sum_j choose(a, j) B(j + 1, a - j + b + 1) B(c + a - j + 1, b + d + 1),
#'   from pu = puv + (1 - puv) pui and 1 - pu = (1 - puv)(1 - pui).
#' Model 3 (puv = pc): sum_j choose(a, j) B(nc + j + 1, 2a + 2b - j + 1)
#'   B(c + a - j + 1, b + d + 1).
#' @keywords internal
#' @noRd
kodikara2020_posterior <- function(yc, yu, model) {
  T <- length(yc)
  nc <- sum(yc)
  tn <- max(which(yc == 1L))
  Nu <- sum(yu)
  cu <- c(0, cumsum(yu))
  k <- seq(tn + 1, T + 1)
  k1 <- pmin(k, T + 1) - 1
  a <- cu[k1 + 1]
  b <- k1 - nc - a
  cc <- Nu - a
  d <- T - k1 - cc
  ll <- vapply(seq_along(k), function(i) {
    if (model == 1) return(lbeta(nc + 1, k1[i] - nc + 1))
    j <- seq(0, a[i])
    lc <- lchoose(a[i], j)
    tail_part <- lbeta(cc[i] + a[i] - j + 1, b[i] + d[i] + 1)
    if (model == 2) {
      lbeta(nc + 1, a[i] + b[i] + 1) +
        kodikara2020_lse(lc + lbeta(j + 1, a[i] - j + b[i] + 1) + tail_part)
    } else {
      kodikara2020_lse(lc + lbeta(nc + j + 1, 2 * a[i] + 2 * b[i] - j + 1) +
                         tail_part)
    }
  }, numeric(1))
  # marginal prior of tau_E: int theta (1 - theta)^(k - 1) d theta
  lp <- ifelse(k <= T, -log(k) - log(k + 1), -log(T + 1))
  lw <- ll + lp
  w <- exp(lw - max(lw))
  list(k = k, w = w / sum(w), pc_a = rep(nc + 1, length(k)),
       pc_b = if (model == 1) k1 - nc + 1 else a + b + 1)
}

#' Median and HDI of tau_E (as period indices) from finite weights + tail
#'
#' For k > T, P(tau_E = k) = w_tail (T + 1) / (k (k + 1)), so the CDF is
#' F(k) = 1 - w_tail (T + 1) / (k + 1).
#' @keywords internal
#' @noRd
kodikara2020_summary <- function(w, k, T, level) {
  nf <- length(k) - 1L
  tail <- w[nf + 1L]
  kf <- k[seq_len(nf)]
  cw <- cumsum(w[seq_len(nf)])
  # smallest index with F >= p
  quant <- function(p) {
    i <- which(cw >= p - 1e-12)
    if (length(i)) return(kf[i[1]])
    if (1 - p <= 1e-12 || tail <= 0) return(Inf)
    max(T + 1, ceiling(tail * (T + 1) / (1 - p) - 1 - 1e-9))
  }
  med <- quant(0.5)
  best <- c(NA, NA)
  width <- Inf
  for (i in seq_len(nf + 1L)) {
    base <- if (i == 1L) 0 else cw[i - 1L]
    start <- if (i <= nf) kf[i] else T + 1
    up <- quant(base + level)
    if (!is.finite(up)) next
    if (up - start < width) {
      width <- up - start
      best <- c(start, up)
    }
  }
  c(median = med, lower = best[1], upper = best[2])
}

#' Posterior median and HDI of pc, a mixture of Beta(pc_a, pc_b)
#' @keywords internal
#' @noRd
kodikara2020_pc <- function(w, pc_a, pc_b, level) {
  x <- seq(0, 1, length.out = 5001)
  dens <- numeric(length(x))
  for (i in seq_along(w)) {
    if (w[i] > 1e-14) dens <- dens + w[i] * stats::dbeta(x, pc_a[i], pc_b[i])
  }
  mass <- dens / sum(dens)
  ord <- order(mass, decreasing = TRUE)
  inside <- ord[seq_len(which(cumsum(mass[ord]) >= level)[1])]
  med <- stats::uniroot(function(p) {
    sum(w * stats::pbeta(p, pc_a, pc_b)) - 0.5
  }, c(0, 1), tol = 1e-10)$root
  c(median = med, lower = x[min(inside)], upper = x[max(inside)])
}
