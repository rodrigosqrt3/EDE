#' Test and confidence interval for a common end-point of several taxa
#' (Solow 1996)
#'
#' Likelihood ratio test of whether several taxa share the same upper
#' end-point, such as the extinction of all of them in one event, and the
#' confidence interval for the common end-point (Solow, 1996). Each taxon's
#' finds are taken to be a stationary Poisson process, so that, given their
#' number, they are independent uniform draws over the unknown range of the
#' taxon. The test improves on the conservative test of Strauss & Sadler
#' (1989) because its significance level is known exactly.
#'
#' For taxon `j` with `n[j]` finds, lowest find `L[j]`, highest find `H[j]`,
#' and `Z` the highest find of all, let
#' `Q[j] = ((H[j] - L[j]) / (Z - L[j]))^(n[j] - 1)`. The likelihood ratio
#' statistic, conditional on the lowest find of each taxon, is the product of
#' the `Q[j]`, and under the null hypothesis of a common upper end-point
#' `-2 log(product)` has a chi-square distribution with `2 (k - 1)` degrees of
#' freedom (Hogg, 1956), for `k` taxa; with two taxa the statistic is itself
#' the p-value. The confidence interval for the common end-point `U`, assuming
#' it holds, is `(Z, U_u)`, where `U_u` is the value at which `-2 log` of the
#' product of `((H[j] - L[j]) / (U - L[j]))^(n[j] - 1)` equals the upper `alpha`
#' quantile of the chi-square distribution with `2 k` degrees of freedom; it is
#' found numerically. For one taxon it is the interval of [strauss1989()].
#'
#' The upper end-point is the one at the largest time. For a common lower
#' end-point (for instance a common time of first appearance) use
#' `end = "lower"`, which, as the paper recommends, multiplies the locations by
#' minus one. The paper's assumption of independent stationary Poisson
#' processes is crucial and should be checked (see [vogel2009()]); the
#' interval is conditional on a common end-point and becomes too narrow if one
#' taxon ended earlier than the others.
#'
#' @param taxa The taxa, as a data frame with one row per taxon and columns
#'   `n` (number of finds, at least 2), `lowest` and `highest` (the locations
#'   of the lowest and the highest find), or a list of [sighting_data] objects,
#'   one per taxon, in which each distinct time with a positive count is one
#'   find (counts are ignored, as in [strauss1989()]).
#' @param alpha Significance level, in (0, 1), of the confidence interval.
#' @param end `"upper"` (the default) or `"lower"`: which end-point is common.
#'
#' @return An [ede_estimate] object (class `ede_solow1996`) with `estimate =
#'   NA` and the confidence interval for the common end-point in `lower` and
#'   `upper` (for `end = "upper"`, `lower` is the highest find), plus
#'   `statistic` (`-2 log` of the likelihood ratio statistic), `df`, `p_value`
#'   (`NA` for a single taxon) and `taxa` (the data used, in the original
#'   orientation).
#'
#' @references
#' Solow, A. R. (1996). Tests and confidence intervals for a common upper
#'   endpoint in fossil taxa. Paleobiology, 22(3), 406-410.
#'
#' Strauss, D., & Sadler, P. M. (1989). Classical confidence intervals and
#'   Bayesian probability estimates for ends of local taxon ranges.
#'   Mathematical Geology, 21(4), 411-427.
#'
#' @seealso [strauss1989()] for one taxon, [seymour_ammonites].
#'
#' @examples
#' # ten uppermost ammonite taxa of Seymour Island (Solow 1996, Table 1)
#' solow1996(seymour_ammonites)
#' # two of them, with a 90% interval
#' solow1996(seymour_ammonites[c(3, 8), ], alpha = 0.1)
#'
#' @export
solow1996 <- function(taxa, alpha = 0.05, end = c("upper", "lower")) {
  check_alpha(alpha)
  end <- match.arg(end)
  d <- s96_taxa(taxa)
  orig <- d
  if (end == "lower") {
    d <- data.frame(n = orig$n, lowest = -orig$highest,
                    highest = -orig$lowest)
  }
  k <- nrow(d)
  z <- max(d$highest)
  # -2 log of the product of ((H - L) / (U - L))^(n - 1) over the taxa
  q <- function(u) {
    -2 * sum((d$n - 1) * log((d$highest - d$lowest) / (u - d$lowest)))
  }
  stat <- q(z)
  df <- 2 * (k - 1)
  p <- if (k > 1L) stats::pchisq(stat, df, lower.tail = FALSE) else NA_real_
  crit <- stats::qchisq(1 - alpha, 2 * k)
  if (stat >= crit) {
    warning("a common end-point is rejected at this alpha even at the ",
            "highest find; the interval is empty and is returned as that ",
            "find.", call. = FALSE)
    u <- z
  } else {
    hi <- z + max(z - min(d$lowest), 1)
    while (q(hi) < crit) hi <- z + 2 * (hi - z)
    u <- stats::uniroot(function(x) q(x) - crit, c(z, hi), tol = 1e-10)$root
  }
  if (end == "upper") {
    lower <- z
    upper <- u
  } else {
    lower <- -u
    upper <- -z
  }
  out <- new_ede_estimate(NA_real_, lower, upper,
                          method = "Solow (1996)", alpha = alpha,
                          interval_type = "one-sided")
  out$statistic <- stat
  out$df <- df
  out$p_value <- p
  out$end <- end
  out$taxa <- orig
  class(out) <- c("ede_solow1996", class(out))
  out
}

#' @export
print.ede_solow1996 <- function(x, ...) {
  cat(sprintf("<%s, common %s end-point>\n", x$method, x$end))
  cat(sprintf("  %d taxa\n", nrow(x$taxa)))
  if (is.na(x$p_value)) {
    cat("  test of a common end-point: not defined for one taxon\n")
  } else {
    cat(sprintf("  test of a common end-point: -2 log L = %s on %d df, p = %s\n",
                format(round(x$statistic, 2)), as.integer(x$df),
                format(signif(x$p_value, 3))))
  }
  cat(sprintf("  %g%% confidence interval for the end-point: [%s, %s]\n",
              100 * (1 - x$alpha), format(round(x$lower, 2)),
              format(round(x$upper, 2))))
  invisible(x)
}

#' Summary of the finds of each taxon
#' @keywords internal
#' @noRd
s96_taxa <- function(taxa) {
  if (inherits(taxa, "sighting_data")) {
    stop("`taxa` needs one entry per taxon: a data frame of `n`, `lowest` ",
         "and `highest`, or a list of sighting_data objects.", call. = FALSE)
  }
  if (is.data.frame(taxa)) {
    if (!all(c("n", "lowest", "highest") %in% names(taxa))) {
      stop("`taxa` must have columns `n`, `lowest` and `highest`.",
           call. = FALSE)
    }
    d <- data.frame(n = taxa$n, lowest = taxa$lowest, highest = taxa$highest)
  } else if (is.list(taxa) && length(taxa) > 0L &&
             all(vapply(taxa, inherits, logical(1), "sighting_data"))) {
    rows <- lapply(taxa, function(sd) {
      t <- unique(sd$time[sd$count > 0])
      if (!length(t)) stop("a taxon has no sightings with count > 0.",
                           call. = FALSE)
      data.frame(n = length(t), lowest = min(t), highest = max(t))
    })
    d <- do.call(rbind, rows)
  } else {
    stop("`taxa` must be a data frame or a list of sighting_data objects.",
         call. = FALSE)
  }
  if (!is.numeric(d$n) || !is.numeric(d$lowest) || !is.numeric(d$highest) ||
      anyNA(d) || any(!is.finite(as.matrix(d)))) {
    stop("`n`, `lowest` and `highest` must be finite numbers.", call. = FALSE)
  }
  if (any(d$n != round(d$n)) || any(d$n < 2)) {
    stop("each taxon needs at least 2 finds (`n` a whole number).",
         call. = FALSE)
  }
  if (any(d$highest <= d$lowest)) {
    stop("`highest` must be greater than `lowest` for every taxon.",
         call. = FALSE)
  }
  rownames(d) <- NULL
  d
}
