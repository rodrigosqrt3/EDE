#' Gaussian-resampled, inverse-weighted McInerny estimate (GRIWM)
#' (Bradshaw et al. 2012)
#'
#' Estimate of the extinction time from a series of dated records that allows
#' for dating error, by the GRIWM method of Bradshaw et al. (2012). It extends
#' the sighting-rate model of McInerny et al. (2006), see [mcinerny2006()], in
#' two ways. First, the estimate is computed for every number `k` of the most
#' recent records, from 2 to `n`, and the results are averaged with weights
#' that decrease with the age of the `k`-th record, because the most recent
#' records are taken to be more informative about the sighting rate as
#' extinction is approached. Second, each date is resampled from a normal
#' distribution with its dating error and the weighted estimate is recomputed
#' for each draw.
#'
#' For the `k` most recent records, with `x[1]` the most recent and
#' `s[k] = x[1] - x[k]` the span of the series, the probability of no further
#' record up to a time `x[1] + d` after the last one is `(1 - k / s[k])^d`
#' (equation 3 of the
#' paper, with `k` records over a period `s[k]`). The estimate `theta[k]` for
#' level `alpha` is the first time at which this probability is `alpha`,
#' `x[1] + log(alpha) / log(1 - k / s[k])`. The weighted estimate is
#' `sum(w[k] theta[k]) / sum(w[k])` over `k = 2, ..., n`, with weights
#' `w[k] = 1 / s[k]` (equation 4; the paper expresses them as a ratio to a
#' common constant, which cancels). The time unit is the one in which the
#' sighting probability per unit is meant, as in [mcinerny2006()]: for
#' radiocarbon ages, years, not thousands of years.
#'
#' With dating errors, the estimate is the median of `nsim` draws of the
#' weighted estimate and the interval is formed by the quantiles of the draws,
#' as in the paper. Without errors, the result is the inverse-weighted
#' McInerny estimate of the paper (its Figure 2) and no interval is given. For
#' a draw, `theta[k]` is undefined if the series has `k` or more records in a
#' span of at most `k` time units (the probability above is then not defined);
#' such values are left out of the weighted mean, and a draw with none left is
#' dropped and counted in `failed`.
#'
#' The paper compares the method with others on simulated series and reports
#' that it has the lowest average bias, with wider intervals, when the
#' sighting rate declines towards extinction; the weights are an heuristic of
#' the authors, not a derived estimator.
#'
#' @param sd A [sighting_data] object whose positive-count rows are the dated
#'   records (a count above one repeats the date), or, with `error`, a data
#'   frame with columns `time` and `error`, one row per record.
#' @param error Dating standard deviation of each record: one value, or one per
#'   row of `sd` in the order of its rows; or `NULL` for no dating error. For a
#'   data frame with an `error` column it is taken from the data frame.
#' @param alpha Threshold probability of the sighting-rate model, in (0, 1);
#'   the paper uses 0.05.
#' @param level Coverage of the interval from the resampled estimates, in
#'   (0, 1); the paper uses 0.95.
#' @param nsim Number of resampled series; the paper uses 10,000.
#' @param k Number of most recent records to use (at least 2); by default
#'   all. The paper uses the 25 most recent.
#'
#' @return An [ede_estimate] object with the median weighted estimate in
#'   `estimate` and, if there are dating errors, the interval in `lower` and
#'   `upper` (`interval_type` `"two-sided"`), and the components `theta_k` (a
#'   data frame with `k`, `span` and `theta` for the dates as recorded),
#'   `theta_w` (the weighted estimate for the dates as recorded), `draws` and
#'   `failed`.
#'
#' @references
#' Bradshaw, C. J. A., Cooper, A., Turney, C. S. M., & Brook, B. W. (2012).
#'   Robust estimates of extinction time in the geological record. Quaternary
#'   Science Reviews, 33, 14-19. doi:10.1016/j.quascirev.2011.11.021
#'
#' McInerny, G. J., Roberts, D. L., Davy, A. J., & Cribb, P. J. (2006).
#'   Significance of sighting rate in inferring extinction and threat.
#'   Conservation Biology, 20(2), 562-567.
#'
#' @seealso [mcinerny2006()], [solow2006()] for the likelihood approach to
#'   dated records, [alaska_megafauna].
#'
#' @examples
#' mam <- alaska_megafauna[alaska_megafauna$taxon == "mammoth", ]
#' # years, increasing towards the present
#' dates <- data.frame(time = -mam$age_ka_bp * 1000, error = mam$sd_ka * 1000)
#' set.seed(1)
#' bradshaw2012(dates, nsim = 1000)
#' # without dating error: the inverse-weighted McInerny estimate
#' bradshaw2012(sighting_data(data.frame(time = dates$time, count = 1)))
#'
#' @export
bradshaw2012 <- function(sd, error = NULL, alpha = 0.05, level = 0.95,
                         nsim = 10000, k = NULL) {
  check_alpha(alpha)
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      level <= 0 || level >= 1) {
    stop("`level` must be in (0, 1).", call. = FALSE)
  }
  if (!is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) || nsim < 1) {
    stop("`nsim` must be a positive whole number.", call. = FALSE)
  }
  nsim <- as.integer(nsim)
  d <- gw_data(sd, error, k)
  x <- d$x
  s <- d$s

  fit <- gw_estimate(x, alpha)
  theta_k <- data.frame(k = seq_along(x), span = x[1L] - x,
                        theta = fit$theta)[-1L, , drop = FALSE]
  has_error <- any(s > 0)

  if (!has_error) {
    if (is.na(fit$tw)) {
      stop("the estimate is undefined: no number of records `k` has fewer ",
           "records than time units in its span.", call. = FALSE)
    }
    out <- new_ede_estimate(fit$tw, method = "Bradshaw et al. (2012), no dating error",
                            alpha = 1 - level)
    out$theta_k <- theta_k
    out$theta_w <- fit$tw
    return(out)
  }

  draws <- vapply(seq_len(nsim), function(i) {
    gw_estimate(stats::rnorm(length(x), x, s), alpha)$tw
  }, numeric(1))
  ok <- !is.na(draws)
  failed <- sum(!ok)
  if (!any(ok)) {
    stop("the estimate is undefined for every resampled series.",
         call. = FALSE)
  }
  if (failed > 0L) {
    warning(failed, " of ", nsim, " resampled series gave no estimate and ",
            "were dropped.", call. = FALSE)
  }
  q <- stats::quantile(draws[ok], c((1 - level) / 2, 1 - (1 - level) / 2),
                       names = FALSE)
  out <- new_ede_estimate(stats::median(draws[ok]), q[1L], q[2L],
                          method = "Bradshaw et al. (2012)",
                          alpha = 1 - level)
  out$theta_k <- theta_k
  out$theta_w <- fit$tw
  out$draws <- draws[ok]
  out$failed <- failed
  out
}

#' Records, most recent first, with their dating errors
#' @keywords internal
#' @noRd
gw_data <- function(sd, error, k) {
  if (is.null(error) && inherits(sd, "sighting_data")) {
    pos <- sd$count > 0
    if (!any(pos)) stop("no sightings with count > 0.", call. = FALSE)
    x <- rep(sd$time[pos], sd$count[pos])
    s <- rep(0, length(x))
    o <- order(x, decreasing = TRUE)
    x <- x[o]
    s <- s[o]
  } else {
    d <- sr_data(sd, error, k = NULL, kmin = 2L)
    x <- d$x
    s <- d$s
  }
  n <- length(x)
  if (n < 2L) stop("at least 2 dated records are needed.", call. = FALSE)
  if (is.null(k)) k <- n
  if (length(k) != 1L || !is.numeric(k) || is.na(k) || k != as.integer(k) ||
      k < 2L || k > n) {
    stop("`k` must be an integer between 2 and the number of records (", n,
         ").", call. = FALSE)
  }
  list(x = x[seq_len(k)], s = s[seq_len(k)])
}

#' McInerny estimates for each k and their inverse-span weighted mean
#'
#' `x` is sorted into decreasing order here, so a resampled series is allowed.
#' @keywords internal
#' @noRd
gw_estimate <- function(x, alpha) {
  x <- sort(x, decreasing = TRUE)
  n <- length(x)
  k <- seq_len(n)
  span <- x[1L] - x
  theta <- rep(NA_real_, n)
  idx <- k[-1L]
  ok <- span[idx] > idx # k / span < 1
  theta[idx[ok]] <- x[1L] + log(alpha) / log(1 - idx[ok] / span[idx[ok]])
  w <- 1 / span
  use <- !is.na(theta)
  tw <- if (any(use)) sum(w[use] * theta[use]) / sum(w[use]) else NA_real_
  list(theta = theta, tw = tw)
}
