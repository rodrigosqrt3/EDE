#' Choose the number of recent sightings for extreme-value estimators
#'
#' Extreme-value estimators such as [ole()] and [hall_wang1999()] use only the
#' `k` most recent sightings, because only the extreme order statistics of the
#' sighting times follow the Weibull extreme-value model. Too large a `k`
#' biases the estimate; too small a `k` makes it imprecise. `ede_select_k()`
#' chooses `k` automatically, adapting to this problem the threshold selection
#' of Bader, Yan & Zhang (2018): goodness-of-fit tests of the Weibull model
#' are run for `k = kmax, kmax - 1, ..., kmin`, from the most to the least
#' data, and the ForwardStop rule of G'Sell et al. (2016) controls the false
#' discovery rate of the ordered sequence of tests. The chosen `k` is the
#' largest one that is not rejected.
#'
#' The test for a given `k` is a parametric bootstrap. Its statistic is the
#' minimum over candidate end-points of Greenwood's statistic of the
#' log-spacings (the criterion of [hall_wang1999()]), which is large when the
#' `k` latest sightings do not follow the model. The null distribution is
#' simulated from the Weibull model fitted to the same `k` sightings (extremal
#' exponent estimated by the minimum-distance method), so each p-value is
#' approximate: it ignores the uncertainty of the fitted exponent.
#'
#' This is a package construct that applies the ordered-testing idea of Bader
#' et al. (2018), developed for generalised Pareto thresholds, to the
#' sighting-record model; it is not a published method for sighting records.
#' The procedure is random: call [set.seed()] first for reproducible results.
#'
#' @inheritParams hall_wang1999
#' @param kmin,kmax Smallest and largest `k` examined. `kmin` is at least 4;
#'   `kmax` defaults to the number of sightings, at most 30.
#' @param alpha False discovery rate of the ForwardStop rule, in (0, 1).
#' @param nsim Number of bootstrap draws per value of `k`.
#'
#' In simulations from the Weibull model the single tests reject at close to
#' their nominal level and the procedure keeps the largest `k` in about nine
#' cases out of ten; against a record whose older sightings come from a much
#' more diffuse regime its power is moderate. A `k` for which the
#' minimum-distance estimator of [hall_wang1999()] does not exist can still be
#' selected, because a record with widely spread older sightings is compatible
#' with the Gumbel limit of the model; the `estimable` column of the table
#' flags these cases.
#'
#' @return An object of class `ede_select_k`: a list with `k` (the selected
#'   value, `NA` if every `k` is rejected), `table` (a data frame with columns
#'   `k`, `statistic`, `shape`, `p_value`, `estimable` and `rejected`),
#'   `alpha` and `method`.
#'
#' @references
#' Bader, B., Yan, J., & Zhang, X. (2018). Automated threshold selection for
#'   extreme value analysis via ordered goodness-of-fit tests with adjustment
#'   for false discovery rate. The Annals of Applied Statistics, 12(1),
#'   310-329. doi:10.1214/17-AOAS1092
#'
#' G'Sell, M. G., Wager, S., Chouldechova, A., & Tibshirani, R. (2016).
#'   Sequential selection procedures and false discovery rate control. Journal
#'   of the Royal Statistical Society B, 78(2), 423-444.
#'
#' @seealso [hall_wang1999()], [ole()]
#'
#' @examples
#' sd <- sighting_data(data.frame(year = c(1880, 1900, 1904, 1909, 1911, 1913,
#'                                          1918, 1920, 1925, 1928, 1930, 1934),
#'                                sightings = 1))
#' set.seed(1)
#' sel <- ede_select_k(sd, nsim = 99)
#' sel
#'
#' @export
ede_select_k <- function(sd, kmin = 4L, kmax = NULL, alpha = 0.1,
                         events = FALSE, nsim = 499L) {
  stopifnot(inherits(sd, "sighting_data"))
  check_alpha(alpha)
  if (!is.numeric(nsim) || length(nsim) != 1L || is.na(nsim) || nsim < 19) {
    stop("`nsim` must be a number of at least 19.", call. = FALSE)
  }
  t <- evt_times(sd, events)
  n <- length(t)
  kmin <- evt_k(kmin, n, kmin = 4L, default = 4L, label = "ede_select_k()")
  if (is.null(kmax)) kmax <- min(n, 30L)
  kmax <- evt_k(kmax, n, kmin = kmin, default = kmax, label = "ede_select_k()")

  ks <- kmax:kmin
  rows <- lapply(ks, function(k) {
    g <- t[1L] - t[seq_len(k)]
    if (g[k] <= 0) {
      return(data.frame(k = k, statistic = NA_real_, shape = NA_real_,
                        p_value = NA_real_, estimable = NA))
    }
    fit <- hw_fit(g)
    stat <- hw_min_stat(g)
    shape <- if (is.na(fit$shape)) {
      # no interior minimum: fit the exponent at the edge of the search grid
      s <- max(hw_grid(g))
      (k - 1L) / sum(seq_len(k - 1L) * diff(log(s + g)))
    } else {
      fit$shape
    }
    sim <- vapply(seq_len(nsim), function(b) {
      gg <- cumsum(stats::rexp(k))^(1 / shape)
      hw_min_stat(gg - gg[1L])
    }, numeric(1))
    data.frame(k = k, statistic = stat, shape = shape,
               p_value = (1 + sum(sim >= stat)) / (nsim + 1),
               estimable = !is.na(fit$s))
  })
  tab <- do.call(rbind, rows)
  if (anyNA(tab$p_value)) {
    stop("the `k` latest sightings do not span more than one time for some ",
         "`k`; use events = FALSE or a larger kmin.", call. = FALSE)
  }

  j <- forward_stop(tab$p_value, alpha)
  tab$rejected <- seq_len(nrow(tab)) <= j
  k_sel <- if (j >= nrow(tab)) NA_integer_ else tab$k[j + 1L]
  if (is.na(k_sel)) {
    warning("the Weibull model is rejected for every `k` examined; ",
            "returning NA.", call. = FALSE)
  }
  structure(list(k = k_sel, table = tab, alpha = alpha,
                 method = "Ordered goodness-of-fit selection of k"),
            class = "ede_select_k")
}

#' ForwardStop rule of G'Sell et al. (2016)
#'
#' Number of leading hypotheses rejected: the largest `j` for which the mean
#' of `-log(1 - p)` over the first `j` ordered p-values is at most `alpha`.
#' @keywords internal
#' @noRd
forward_stop <- function(p, alpha) {
  cum <- cumsum(-log1p(-p)) / seq_along(p)
  ok <- which(cum <= alpha)
  if (length(ok)) max(ok) else 0L
}

#' @export
print.ede_select_k <- function(x, ...) {
  cat(sprintf("<%s>\n", x$method))
  cat(sprintf("  selected k: %s  (ForwardStop, FDR %g)\n",
              if (is.na(x$k)) "none" else format(x$k), x$alpha))
  tab <- x$table
  tab$statistic <- round(tab$statistic, 4)
  tab$shape <- round(tab$shape, 2)
  tab$p_value <- round(tab$p_value, 3)
  print(tab, row.names = FALSE)
  invisible(x)
}

#' @export
plot.ede_select_k <- function(x, ...) {
  tab <- x$table
  graphics::plot(tab$k, tab$p_value, type = "b", ylim = c(0, 1),
                 xlab = "k (most recent sightings)",
                 ylab = "goodness-of-fit p-value", pch = ifelse(tab$rejected, 4, 16),
                 ...)
  graphics::abline(h = x$alpha, lty = 2, col = "grey50")
  if (!is.na(x$k)) graphics::abline(v = x$k, lty = 3, col = "firebrick")
  invisible(x)
}
