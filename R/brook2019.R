#' Reliability-weighted resampling of a sighting record (Brook et al., 2019)
#'
#' Applies any extinction-date estimator to many random sub-records drawn
#' from a sighting record of mixed reliability, and summarises the resulting
#' empirical distribution (the "BBJ" weighted-resampling approach of Brook,
#' Buettel & Jaric, 2019).
#'
#' The algorithm follows the paper and the authors' R code (Data S1, Zenodo
#' doi:10.5281/zenodo.3247189, `bbj.2018()`):
#' 1. Sightings in the same time unit are combined into one record with
#'    reliability `p[j] = 1 - prod(1 - p[ij])` (the unnumbered equation in
#'    the Methods of Brook et al., 2019).
#' 2. Each record is kept independently with probability equal to its
#'    reliability (resampling "without replacement", using reliability as an
#'    inclusion probability). Records with reliability 1 are therefore always
#'    kept; there is no other rule that forces the most recent certain
#'    sighting to be retained.
#' 3. `method` is applied to each sub-record and the chosen `statistic` is
#'    stored. Sub-records with fewer than `min_sightings` records (the
#'    authors' code requires at least 3), and sub-records for which `method`
#'    fails or returns `NA`, are discarded; their number and proportion are
#'    reported.
#' 4. The kept values form the empirical distribution of the statistic.
#'
#' Summaries: `estimate` is the median of the distribution, the point
#' estimate named in the paper's Methods and in its Fig. 1 (the authors' code
#' instead reports the mean, which is returned as `mean`). `upper` is the
#' `1 - alpha` quantile of the distribution, the "upper 95% confidence bound"
#' of the paper (its text calls it the "2.5th percentile", apparently in
#' error; the authors' code takes the last year at which the proportion of
#' values not yet reached exceeds `alpha`, i.e. the `1 - alpha` quantile).
#' The extinction years in Table 1(b) of the paper are these upper bounds.
#' If `persist_year` is given, `persistence` is the proportion of values
#' later than `persist_year`, the "probability of persistence" in Table 1.
#' Quantiles use `type = 1` (inverse empirical distribution), so infinite
#' values (e.g. a non-existent upper bound) are handled.
#'
#' Kept records enter the sub-record with count 1. Zero-count rows of `sd`
#' (for example a leading row that fixes the observation origin) are kept in
#' every sub-record; otherwise the earliest kept sighting becomes the origin
#' for estimators that use one, as in the authors' code.
#'
#' Status of the method: Solow (2023) criticised the statistical properties
#' of the extinction date and interval produced this way, noting in
#' particular that the estimated extinction time exceeds the true one when
#' reliabilities are high. Jaric et al. (2023) replied that the procedure is
#' not itself an extinction-date estimator but a pre-processing step for
#' the input data, whose result is passed to any estimator (McCrea et al.,
#' 2024). The summaries returned here are descriptive statistics of that
#' resampling distribution, not confidence intervals with established
#' coverage.
#'
#' @param sd A [sighting_data] object.
#' @param reliability Sighting reliabilities in `[0, 1]`. Either a numeric
#'   vector with one value per row of `sd` or per positive-count row (taken
#'   as the combined reliability of that time unit), or a list with one
#'   numeric vector per positive-count row giving the reliabilities of the
#'   individual sightings in that time unit (its length must equal the
#'   count); these are combined as in step 1. Values for zero-count rows are
#'   ignored.
#' @param method An estimator function whose first argument is a
#'   [sighting_data] object, such as [ole()], [strauss1989()],
#'   [robson1964()] or [solow1993()]. It may return an [ede_estimate] object
#'   or a single number. The "S93" estimator of the authors' code,
#'   `t[n] + (t[n] - t[1]) / (n - 1)` with the first sighting as origin,
#'   is the point estimate of [strauss1989()] (the default).
#' @param n_resamples Number of random sub-records (the paper suggests
#'   about 1,000 or more; the authors' code uses 10,000).
#' @param statistic Component of the [ede_estimate] returned by `method` to
#'   collect: `"estimate"` (default) or `"upper"`. Ignored if `method`
#'   returns a number.
#' @param alpha Level for the upper quantile, in (0, 1).
#' @param persist_year Optional time at which to report the proportion of
#'   resampled values later than it.
#' @param min_sightings Minimum number of kept records for a sub-record to be
#'   analysed.
#' @param seed Optional random seed. The global random number state is
#'   restored on exit.
#' @param ... Further arguments passed to `method` (e.g. `test_year`).
#'
#' @return An [ede_estimate] object (`interval_type = "one-sided"`) with
#'   `estimate` (median), `upper` (`1 - alpha` quantile) and additional
#'   components `mean`, `values` (the kept resampled values),
#'   `n_resamples`, `n_failed`, `prop_failed`, `statistic`,
#'   `combined_reliability` and, if `persist_year` is given, `persistence`
#'   and `persist_year`.
#'
#' @references
#' Brook, B. W., Buettel, J. C., & Jaric, I. (2019). A fast re-sampling
#'   method for using reliability ratings of sightings with extinction-date
#'   estimators. Ecology, 100(9), e02787.
#'
#' Jaric, I., Buettel, J. C., & Brook, B. W. (2023). A fast re-sampling
#'   method for using reliability ratings of sightings with extinction-date
#'   estimators: reply. Ecology, 104(9), e4124.
#'
#' Solow, A. R. (2023). A fast re-sampling method for using reliability
#'   ratings of sightings with extinction-date estimators: comment. Ecology,
#'   104(9), e4123.
#'
#' McCrea, R. S., Cheale, T., Campillo-Funollet, E., & Roberts, D. L.
#'   (2024). Inferring species extinction from sighting data. Cambridge
#'   Prisms: Extinction, 2, e19.
#'
#' @export
brook2019 <- function(sd, reliability, method = strauss1989,
                      n_resamples = 1000, statistic = c("estimate", "upper"),
                      alpha = 0.05, persist_year = NULL, min_sightings = 3,
                      seed = NULL, ...) {
  stopifnot(inherits(sd, "sighting_data"))
  statistic <- match.arg(statistic)
  method <- match.fun(method)
  if (alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0, 1).", call. = FALSE)
  if (length(n_resamples) != 1L || !is.numeric(n_resamples) ||
      is.na(n_resamples) || n_resamples < 1 ||
      n_resamples != floor(n_resamples)) {
    stop("`n_resamples` must be a positive whole number.", call. = FALSE)
  }
  if (length(min_sightings) != 1L || !is.numeric(min_sightings) ||
      is.na(min_sightings) || min_sightings < 1) {
    stop("`min_sightings` must be a number >= 1.", call. = FALSE)
  }
  if (!is.null(persist_year) && (length(persist_year) != 1L ||
      !is.numeric(persist_year) || is.na(persist_year))) {
    stop("`persist_year` must be a single number.", call. = FALSE)
  }
  occupied <- which(sd$count > 0)
  if (length(occupied) == 0L) stop("no sightings with count > 0.", call. = FALSE)
  if (missing(reliability)) stop("`reliability` is required.", call. = FALSE)
  p <- brook2019_combine(reliability, sd, occupied)

  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = globalenv())
    on.exit({
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = globalenv())
      } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        rm(".Random.seed", envir = globalenv())
      }
    }, add = TRUE)
    set.seed(seed)
  }

  zero_times <- sd$time[sd$count == 0]
  zero_rows <- data.frame(time = zero_times, count = rep(0, length(zero_times)))
  occ_times <- sd$time[occupied]
  one_value <- function() {
    keep <- stats::runif(length(p)) < p
    if (sum(keep) < min_sightings) return(NA_real_)
    sub <- rbind(zero_rows, data.frame(time = occ_times[keep], count = 1))
    sub <- sighting_data(sub)
    res <- tryCatch(suppressWarnings(method(sub, ...)),
                    error = function(e) NULL)
    if (is.null(res)) return(NA_real_)
    val <- if (inherits(res, "ede_estimate")) res[[statistic]] else res
    if (!is.numeric(val) || length(val) != 1L) return(NA_real_)
    as.numeric(val)
  }
  values <- vapply(seq_len(n_resamples), function(i) one_value(), numeric(1))

  failed <- is.na(values)
  values <- values[!failed]
  if (length(values) == 0L) {
    stop("`method` returned no usable value for any resampled record.",
         call. = FALSE)
  }
  q <- function(prob) unname(stats::quantile(values, prob, type = 1))

  inner <- tryCatch({
    full <- sd
    full$count[occupied] <- 1
    m <- suppressWarnings(method(full, ...))
    if (inherits(m, "ede_estimate")) m$method else NULL
  }, error = function(e) NULL)
  label <- if (is.null(inner)) {
    "Brook et al. (2019) resampling"
  } else {
    sprintf("Brook et al. (2019) resampling: %s", inner)
  }

  result <- new_ede_estimate(q(0.5), lower = NA_real_, upper = q(1 - alpha),
                             method = label, alpha = alpha,
                             interval_type = "one-sided")
  result$mean <- mean(values)
  result$values <- values
  result$n_resamples <- n_resamples
  result$n_failed <- sum(failed)
  result$prop_failed <- mean(failed)
  result$statistic <- statistic
  result$combined_reliability <- p
  if (!is.null(persist_year)) {
    result$persist_year <- persist_year
    result$persistence <- mean(values > persist_year)
  }
  if (result$prop_failed > 0) {
    message(sprintf(
      "%d of %d resampled records (%.1f%%) were discarded (too few sightings or `method` failed).",
      result$n_failed, n_resamples, 100 * result$prop_failed
    ))
  }
  result
}

#' Validate reliabilities and combine sightings within a time unit
#' @keywords internal
#' @noRd
brook2019_combine <- function(reliability, sd, occupied) {
  check <- function(x) {
    if (!is.numeric(x) || anyNA(x) || any(!is.finite(x)) || any(x < 0 | x > 1)) {
      stop("`reliability` must contain finite values in [0, 1].", call. = FALSE)
    }
  }
  if (is.list(reliability)) {
    if (length(reliability) != length(occupied)) {
      stop("a `reliability` list needs one element per positive-count row.",
           call. = FALSE)
    }
    lapply(reliability, check)
    lens <- lengths(reliability)
    if (any(lens != sd$count[occupied])) {
      stop("each `reliability` list element must have one value per sighting (its length must equal the count).",
           call. = FALSE)
    }
    return(vapply(reliability, function(x) 1 - prod(1 - x), numeric(1)))
  }
  check(reliability)
  if (length(reliability) == nrow(sd)) {
    reliability[occupied]
  } else if (length(reliability) == length(occupied)) {
    reliability
  } else {
    stop("`reliability` must have one value per row or per positive-count row.",
         call. = FALSE)
  }
}
