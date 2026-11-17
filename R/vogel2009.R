#' Probability plot correlation test of constant sighting rate
#'
#' Tests whether the sighting times of a record are compatible with a constant
#' sighting rate before extinction, that is, with times uniformly distributed
#' on an interval ending at the extinction time (the assumption of
#' [solow1993()] and [solow2005()]), with the probability plot correlation
#' coefficient (PPCC) test of Vogel et al. (2009), who apply it to the sighting
#' records of North American and Hawaiian birds.
#'
#' With the `n` sighting times in increasing order, `t[1] <= ... <= t[n]`, the
#' statistic is the product-moment correlation `r` between `t[i]` and the
#' plotting positions `i / (n + 1)` (equation 14 of the paper). Low values of
#' `r` are evidence against the uniform model. `r` is invariant to the origin
#' and the unit of time, and to the unknown end of the interval, so its null
#' distribution depends only on `n`. The p-value is the probability, under the
#' uniform model, of a correlation at most `r`; it is computed by simulation
#' (`nsim` draws, with a fixed internal seed, so the result is reproducible and
#' the global random number stream is left untouched) in place of the tables of
#' Vogel & Kroll (1989).
#'
#' As Vogel et al. stress, a single short record gives the test little power,
#' so failing to reject is not evidence for the uniform model. For a set of
#' records, [ede_field_test()] combines the p-values.
#'
#' By default each distinct sighting time counts once, as in
#' [hall_wang1999()]; set `events = TRUE` to count every sighting event, in
#' which case ties are present and the model is only approximate.
#'
#' @param sd A [sighting_data] object.
#' @param events If `FALSE` (default), count each distinct sighting time once;
#'   if `TRUE`, count every sighting event.
#' @param nsim Number of draws used for the null distribution.
#'
#' @return An object of class `htest`, with `statistic` (`r`), `p.value`,
#'   `parameter` (the number of sighting times `n`) and `alternative`
#'   (`"less"`: the correlation is smaller than under the uniform model).
#'
#' @references
#' Vogel, R. M., Hosking, J. R. M., Elphick, C. S., Roberts, D. L., & Reed, J.
#'   M. (2009). Goodness of fit of probability distributions for sightings as
#'   species approach extinction. Bulletin of Mathematical Biology, 71(3),
#'   701-719. doi:10.1007/s11538-008-9377-3
#'
#' Vogel, R. M., & Kroll, C. N. (1989). Low-flow frequency analysis using
#'   probability-plot correlation coefficients. Journal of Water Resources
#'   Planning and Management, 115(3), 338-357.
#'
#' @seealso [ede_field_test()] to combine tests over several records;
#'   [solow1993()] and [solow2005()], which assume a constant rate.
#'
#' @examples
#' # Eskimo curlew, the example of Vogel et al. (2009, Figure 2)
#' yrs <- c(1946, 1947, 1948, 1950, 1955, 1956, 1959, 1960, 1961, 1962, 1963,
#'          1964, 1968, 1970, 1972, 1973, 1974, 1976, 1977, 1980, 1981, 1982,
#'          1983, 1985)
#' vogel2009(sighting_data(data.frame(year = yrs, sightings = 1)))
#'
#' @export
vogel2009 <- function(sd, events = FALSE, nsim = 1e5) {
  stopifnot(inherits(sd, "sighting_data"))
  check_nsim(nsim)
  t <- sort(evt_times(sd, events))
  n <- length(t)
  if (n < 4L) {
    stop("the PPCC test needs at least 4 sighting times; the record has ", n,
         ".", call. = FALSE)
  }
  if (!(t[n] > t[1L])) {
    stop("the sighting times must not all be equal.", call. = FALSE)
  }
  r <- ppcc_uniform(t)
  p <- ppcc_uniform_p(r, n, nsim)
  structure(
    list(statistic = c(r = r), parameter = c(n = n), p.value = p,
         alternative = "less",
         method = "PPCC test of uniform sighting times (Vogel et al. 2009)",
         data.name = deparse(substitute(sd))),
    class = "htest"
  )
}

#' Combine the PPCC tests of several records (field significance)
#'
#' Combines the individual p-values of tests applied to several records, for
#' example the [vogel2009()] test of a constant sighting rate on many species,
#' with the four field-significance procedures of Vogel et al. (2009, Section
#' 6): a Bonferroni-type test, the false discovery rate procedure of Benjamini
#' & Hochberg (1995), a binomial test of the number of rejections, and a PPCC
#' test of the uniformity of the p-values themselves. With a single short
#' record each test has little power, which is the reason for combining many.
#' The procedures assume the records, hence the p-values, are independent.
#'
#' @param p Numeric vector of p-values, one per record.
#' @param alpha Significance level of each individual test (used for the number
#'   of rejections in the binomial test and the Bonferroni threshold, with the
#'   field level `alpha`), in (0, 1).
#' @param nsim Number of draws used for the null distribution of the PPCC of
#'   the p-values.
#'
#' @return A data frame with one row per procedure and columns `test`,
#'   `statistic`, `rejections` (the number of individual tests rejected, where
#'   defined) and `p_value` (the field significance level). For the Bonferroni
#'   and false discovery rate procedures `p_value` is `NA`: they report the
#'   records rejected at field level `alpha`.
#'
#' @references
#' Vogel, R. M., Hosking, J. R. M., Elphick, C. S., Roberts, D. L., & Reed, J.
#'   M. (2009). Goodness of fit of probability distributions for sightings as
#'   species approach extinction. Bulletin of Mathematical Biology, 71(3),
#'   701-719. doi:10.1007/s11538-008-9377-3
#'
#' Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate:
#'   a practical and powerful approach to multiple testing. Journal of the
#'   Royal Statistical Society B, 57(1), 289-300.
#'
#' @seealso [vogel2009()]
#'
#' @examples
#' p <- c(0.92, 0.65, 0.08, 0.90, 0.43, 0.017, 0.001, 0.015, 0.74, 0.99)
#' ede_field_test(p)
#'
#' @export
ede_field_test <- function(p, alpha = 0.05, nsim = 1e5) {
  if (!is.numeric(p) || anyNA(p) || any(p < 0 | p > 1) || length(p) < 4L) {
    stop("`p` must be at least 4 p-values in [0, 1].", call. = FALSE)
  }
  check_alpha(alpha)
  check_nsim(nsim)
  m <- length(p)
  ps <- sort(p)
  bonf <- sum(p < alpha / m)
  bh_idx <- which(ps <= seq_len(m) * alpha / m)
  bh <- if (length(bh_idx)) max(bh_idx) else 0L
  x <- sum(p <= alpha)
  p_bin <- stats::pbinom(x - 1L, m, alpha, lower.tail = FALSE)
  r <- ppcc_uniform(ps)
  p_ppcc <- ppcc_uniform_p(r, m, nsim)
  data.frame(
    test = c("Bonferroni-type", "False discovery rate (Benjamini-Hochberg)",
             "Binomial", "PPCC of the p-values"),
    statistic = c(alpha / m, NA_real_, x, r),
    rejections = c(bonf, bh, x, NA_integer_),
    p_value = c(NA_real_, NA_real_, p_bin, p_ppcc),
    stringsAsFactors = FALSE
  )
}

#' Correlation between the ordered values and their plotting positions
#' @keywords internal
#' @noRd
ppcc_uniform <- function(x) {
  x <- sort(x)
  n <- length(x)
  stats::cor(x, seq_len(n) / (n + 1))
}

#' P(R <= r) under uniformity, by simulation
#'
#' The correlation is invariant to scale, so the order statistics of a uniform
#' sample are replaced by the cumulative sums of exponentials. The global
#' random number stream is restored.
#' @keywords internal
#' @noRd
ppcc_cache <- new.env(parent = emptyenv())

ppcc_uniform_null <- function(n, nsim) {
  key <- paste(n, nsim)
  if (!is.null(ppcc_cache[[key]])) return(ppcc_cache[[key]])
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had_seed) old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)
  set.seed(20261113L)
  idx <- seq_len(n) - (n + 1) / 2
  r <- numeric(nsim)
  chunk <- 10000L
  start <- 1L
  while (start <= nsim) {
    b <- min(chunk, nsim - start + 1L)
    s <- matrix(stats::rexp(b * n), b, n)
    for (j in seq_len(n)[-1L]) s[, j] <- s[, j] + s[, j - 1L]
    sc <- s - rowMeans(s)
    r[start:(start + b - 1L)] <- as.vector(sc %*% idx) /
      sqrt(rowSums(sc^2) * sum(idx^2))
    start <- start + b
  }
  out <- sort(r)
  ppcc_cache[[key]] <- out
  out
}

ppcc_uniform_p <- function(r, n, nsim) {
  null <- ppcc_uniform_null(n, nsim)
  findInterval(r, null, left.open = FALSE) / length(null)
}
