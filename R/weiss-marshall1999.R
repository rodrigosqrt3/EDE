#' Weiss & Marshall (1999) Bayesian end point for discretely sampled sections
#'
#' Posterior distribution of the true end point `S` of a stratigraphic range
#' (or of a sighting record) when the section is sampled at discrete points
#' `t[1] < ... < t[J]`, each recorded as presence or absence. While the taxon
#' is extant (`t[j] <= S`) it is found in a sample with probability `p`;
#' after `S` it is never found (equations 1-2 of Weiss & Marshall, 1999). With
#' `a[m]` presences up to the last presence `t[m]` and `b(S)` absences up to
#' `S`, the likelihood is `p^a[m] (1 - p)^b(S)` for `S >= t[m]` (equation 3).
#' The prior on `S` is exponential with mean longevity `l` (equation 4) and the
#' prior on `p` is `Beta(c, d)`. Integrating `p` out gives the step-function
#' posterior of equation 6: on `[t[j], t[j+1])` (with `t[J+1] = Inf`) its
#' mass is `g[j]` of equation 8,
#' `B(a[m] + c, b[j] + d) (exp(-t[j] / l) - exp(-t[j+1] / l))`, with
#' `b[j] = j - a[m]`. Quantiles, moments and the posterior of `p` (a mixture
#' of `Beta(a[m] + c, b[j] + d)` densities, equations 7-8) follow in closed
#' form (Appendix).
#'
#' Each row of `sd` is one sampled horizon: a positive count is a presence and
#' a zero count an absence. Every sampled horizon must be supplied, including
#' the barren samples above the last presence, which carry the information
#' about the end point; unsampled positions must not be given as rows. The
#' earliest row is `t[1]`, where the paper assumes the taxon is already extant
#' (it need not be found there). Counts above one are treated as presence;
#' the model uses presence/absence only.
#'
#' Times may be stratigraphic heights or ages, but the end point estimated is
#' always the one at the largest `time`. To estimate a last appearance from
#' ages in Ma (which decrease upwards), code `time = -age`; to estimate a
#' first appearance from heights, code `time = -height`. `mean_longevity` and
#' the results are on the coded scale.
#'
#' The reported estimate is the posterior median and `upper` is the
#' `1 - alpha` posterior quantile, so `[t[m], upper]` is the minimum-length
#' `1 - alpha` posterior interval the paper reports (`S.95`). The paper
#' tabulates these as distances beyond the last presence (residual lifetime);
#' here they are on the original time scale. Choosing `c` and `d` is part of
#' the analysis: the paper pools the records of similar taxa, see
#' [weiss_marshall1999_prior()].
#'
#' @param sd A [sighting_data] object, one row per sampled horizon (count 0 =
#'   absence).
#' @param prior_p Numeric vector `c(c, d)` of the `Beta(c, d)` prior on the
#'   detection probability `p`; both must be positive.
#' @param mean_longevity Mean `l` of the exponential prior on the time to the
#'   end point (equation 4), in the units of `time`; positive and finite.
#'   The paper uses 800 m (a 21 Myr mean foraminiferal longevity converted
#'   to section thickness).
#' @param alpha Posterior tail probability, in (0, 1), for `upper`.
#' @param test_year Optional time after the last presence; if given, the
#'   posterior probability `P(S > test_year | Y)` that the taxon was still
#'   extant at that time is returned as `persistence_prob`.
#'
#' @return An [ede_estimate] object with `estimate` (posterior median of `S`),
#'   `lower` (last presence `t[m]`), `upper` (`1 - alpha` posterior quantile)
#'   and `interval_type = "one-sided"`, plus components `posterior_mean` and
#'   `posterior_sd` of `S`, `p_mean` and `p_sd` (posterior mean and SD of
#'   `p`), `prior_p`, `mean_longevity` and, if requested, `persistence_prob`.
#'
#' @references
#' Weiss, R. E., & Marshall, C. R. (1999). The uncertainty in the true end
#'   point of a fossil's stratigraphic range when stratigraphic sections are
#'   sampled discretely. Mathematical Geology, 31(4), 435-453.
#'
#' Marshall, C. R. (2010). Using confidence intervals to quantify the
#'   uncertainty in the end-points of stratigraphic ranges. In J. Alroy & G.
#'   Hunt (Eds.), Quantitative Methods in Paleobiology (pp. 291-316). The
#'   Paleontological Society Papers, 16.
#'
#' @export
weiss_marshall1999 <- function(sd, prior_p, mean_longevity, alpha = 0.05,
                               test_year = NULL) {
  stopifnot(inherits(sd, "sighting_data"))
  if (alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0, 1).", call. = FALSE)
  if (missing(prior_p) || !is.numeric(prior_p) || length(prior_p) != 2L ||
      anyNA(prior_p) || any(!is.finite(prior_p)) || any(prior_p <= 0)) {
    stop("`prior_p` must be c(c, d) with both values positive and finite.",
         call. = FALSE)
  }
  if (missing(mean_longevity) || !is.numeric(mean_longevity) ||
      length(mean_longevity) != 1L || is.na(mean_longevity) ||
      !is.finite(mean_longevity) || mean_longevity <= 0) {
    stop("`mean_longevity` must be a positive finite number.", call. = FALSE)
  }
  post <- weiss_marshall1999_posterior(sd, prior_p, mean_longevity)
  tm <- post$t_m

  if (!is.null(test_year)) check_test_year(test_year, tm)

  result <- new_ede_estimate(
    tm + weiss_marshall1999_quantile(post, 0.5),
    lower = tm,
    upper = tm + weiss_marshall1999_quantile(post, 1 - alpha),
    method = "Weiss & Marshall (1999)",
    alpha = alpha,
    interval_type = "one-sided"
  )
  result$interval_kind <- "credible"
  result$posterior_mean <- tm + post$mean
  result$posterior_sd <- post$sd
  result$p_mean <- post$p_mean
  result$p_sd <- post$p_sd
  result$prior_p <- c(c = prior_p[[1]], d = prior_p[[2]])
  result$mean_longevity <- mean_longevity
  if (!is.null(test_year)) {
    result$persistence_prob <- 1 - weiss_marshall1999_cdf(post, test_year - tm)
  }
  result
}

#' Beta prior for p pooled from related records (Weiss & Marshall, 1999)
#'
#' Builds the `Beta(c, d)` prior on the detection probability used by
#' [weiss_marshall1999()] the way Weiss & Marshall (1999, p. 442-443) do for
#' their Eastbourne foraminifera: starting from `Beta(0, 0)`, it adds the
#' presences (`c`) and absences (`d`) of other, comparable taxa from the first
#' sample up to, but excluding, each taxon's last presence (the last presence
#' is excluded because it is a presence by definition).
#'
#' @param records A list of [sighting_data] objects, one per comparable taxon,
#'   each with one row per sampled horizon (count 0 = absence).
#'
#' @return A named numeric vector `c(c = , d = )` to pass as `prior_p`.
#'
#' @references
#' Weiss, R. E., & Marshall, C. R. (1999). The uncertainty in the true end
#'   point of a fossil's stratigraphic range when stratigraphic sections are
#'   sampled discretely. Mathematical Geology, 31(4), 435-453.
#'
#' @export
weiss_marshall1999_prior <- function(records) {
  if (inherits(records, "sighting_data")) records <- list(records)
  if (!is.list(records) || length(records) == 0L) {
    stop("`records` must be a non-empty list of sighting_data objects.", call. = FALSE)
  }
  counts <- vapply(records, function(r) {
    if (!inherits(r, "sighting_data")) {
      stop("every element of `records` must be a sighting_data object.", call. = FALSE)
    }
    y <- r$count > 0
    if (!any(y)) stop("every record needs at least one presence.", call. = FALSE)
    before <- y[seq_len(max(which(y)) - 1L)]
    c(sum(before), sum(!before))
  }, numeric(2))
  out <- rowSums(counts)
  if (any(out <= 0)) {
    warning("pooled prior has a zero parameter; it is improper and ",
            "weiss_marshall1999() will not accept it.", call. = FALSE)
  }
  c(c = out[[1]], d = out[[2]])
}

#' Posterior pieces for Weiss & Marshall (1999), in residual time u = S - t[m]
#' @keywords internal
#' @noRd
weiss_marshall1999_posterior <- function(sd, prior_p, l) {
  y <- sd$count > 0
  if (!any(y)) stop("no presences (count > 0) in the record.", call. = FALSE)
  m <- max(which(y))
  J <- nrow(sd)
  a_m <- sum(y)
  cc <- prior_p[[1]]
  dd <- prior_p[[2]]
  t_m <- sd$time[m]

  j <- m:J
  b_j <- j - a_m
  u_lo <- sd$time[j] - t_m
  u_hi <- c(sd$time[j[-1L]] - t_m, Inf)
  # equation 8 with the common factor exp(-t[m] / l) removed (memoryless prior)
  log_beta <- lbeta(a_m + cc, b_j + dd)
  prior_mass <- exp(-u_lo / l) - exp(-u_hi / l)
  log_g <- log_beta + log(prior_mass)
  w <- exp(log_g - max(log_g))
  w <- w / sum(w)

  # moments of u: integral of u^r exp(-u / l) / l over each interval
  seg_moment <- function(r) {
    l^r * gamma(r + 1) *
      (stats::pgamma(u_hi / l, r + 1) - stats::pgamma(u_lo / l, r + 1)) / prior_mass
  }
  m1 <- sum(w * seg_moment(1))
  m2 <- sum(w * seg_moment(2))

  aa <- a_m + cc
  bb <- b_j + dd
  p1 <- sum(w * aa / (aa + bb))
  p2 <- sum(w * aa * (aa + 1) / ((aa + bb) * (aa + bb + 1)))

  list(t_m = t_m, u_lo = u_lo, u_hi = u_hi, w = w, l = l,
       mean = m1, sd = sqrt(max(m2 - m1^2, 0)),
       p_mean = p1, p_sd = sqrt(max(p2 - p1^2, 0)))
}

#' Posterior CDF of the residual u = S - t[m]
#' @keywords internal
#' @noRd
weiss_marshall1999_cdf <- function(post, u) {
  within <- pmin(pmax(u, post$u_lo), post$u_hi)
  frac <- (exp(-post$u_lo / post$l) - exp(-within / post$l)) /
    (exp(-post$u_lo / post$l) - exp(-post$u_hi / post$l))
  min(1, sum(post$w * frac))
}

#' Posterior quantile of the residual u = S - t[m] (Appendix)
#' @keywords internal
#' @noRd
weiss_marshall1999_quantile <- function(post, prob) {
  q <- cumsum(post$w)
  k <- which(q >= prob)[1L]
  # rounding can leave sum(w) just below prob; take the top of the last interval
  if (is.na(k)) k <- length(q)
  q_prev <- if (k == 1L) 0 else q[k - 1L]
  frac <- min(1, (prob - q_prev) / post$w[k])
  e_lo <- exp(-post$u_lo[k] / post$l)
  e_hi <- exp(-post$u_hi[k] / post$l)
  -post$l * log(e_lo - frac * (e_lo - e_hi))
}
