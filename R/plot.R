#' Plot a sighting record
#'
#' Draws the number of sightings recorded at each time, with the observation
#' origin and the last sighting marked.
#'
#' @param x A [sighting_data] object.
#' @param ... Further arguments passed to [graphics::plot()].
#'
#' @return `x`, invisibly.
#'
#' @examples
#' plot(sighting_data(dodo))
#'
#' @export
plot.sighting_data <- function(x, ...) {
  pos <- x$count > 0
  args <- list(...)
  defaults <- list(xlab = "time", ylab = "sightings", type = "h", lwd = 2,
                   ylim = c(0, max(1, x$count)))
  args <- utils::modifyList(defaults, args)
  do.call(graphics::plot, c(list(x$time[pos], x$count[pos]), args))
  graphics::points(x$time[pos], x$count[pos], pch = 16)
  graphics::abline(v = min(x$time), lty = 3, col = "grey50")
  invisible(x)
}

#' Plot an extinction date estimate
#'
#' Draws the point estimate and its interval on a time axis. An interval
#' without a finite end (`Inf`) is drawn with an arrow. For a persistence test,
#' which has no point estimate, only the upper confidence bound is drawn.
#'
#' @param x An [ede_estimate] object.
#' @param ... Further arguments passed to [graphics::plot()].
#'
#' @return `x`, invisibly.
#'
#' @examples
#' plot(ole(sighting_data(dodo)))
#'
#' @export
plot.ede_estimate <- function(x, ...) {
  vals <- c(x$estimate, x$lower, x$upper)
  fin <- vals[is.finite(vals)]
  if (!length(fin)) stop("nothing to plot: no finite estimate or bound.",
                         call. = FALSE)
  span <- max(diff(range(fin)), 1)
  xlim <- c(min(fin) - 0.1 * span, max(fin) + 0.25 * span)
  args <- utils::modifyList(
    list(xlab = "extinction time", ylab = "", yaxt = "n", xlim = xlim,
         ylim = c(0, 2), main = x$method, type = "n"),
    list(...)
  )
  do.call(graphics::plot, c(list(NA_real_, NA_real_), args))
  lo <- if (is.na(x$lower)) NA_real_ else x$lower
  hi <- if (is.na(x$upper)) NA_real_ else x$upper
  if (!is.na(lo) || !is.na(hi)) {
    a <- if (is.na(lo)) fin[1L] else lo
    b <- if (is.na(hi)) a else if (is.infinite(hi)) xlim[2L] else hi
    graphics::segments(a, 1, b, 1, lwd = 3, col = "steelblue")
    if (isTRUE(is.infinite(hi))) {
      graphics::arrows(xlim[2L] - 0.1 * span, 1, xlim[2L], 1, length = 0.1,
                       col = "steelblue", lwd = 3)
    }
  }
  if (!is.na(x$estimate)) graphics::points(x$estimate, 1, pch = 19, cex = 1.4)
  if (is.na(x$estimate) && is.finite(hi)) {
    graphics::points(hi, 1, pch = 25, bg = "steelblue", cex = 1.4)
  }
  invisible(x)
}

#' Plot the results of [ede_compare()]
#'
#' A forest plot with one row per method: point estimate, interval and, for the
#' persistence tests, the one-sided upper confidence bound (a triangle).
#' Methods that could not be applied are omitted.
#'
#' @param x A result of [ede_compare()].
#' @param ... Further arguments passed to [graphics::plot()].
#'
#' @return `x`, invisibly.
#'
#' @examples
#' plot(ede_compare(sighting_data(monk_seal), test_year = 2050))
#'
#' @export
plot.ede_compare <- function(x, ...) {
  keep <- !(is.na(x$estimate) & is.na(x$lower) & is.na(x$upper))
  d <- x[keep, , drop = FALSE]
  if (!nrow(d)) stop("nothing to plot: no method produced a result.",
                     call. = FALSE)
  vals <- c(d$estimate, d$lower, d$upper)
  fin <- vals[is.finite(vals)]
  span <- max(diff(range(fin)), 1)
  xlim <- c(min(fin) - 0.05 * span, max(fin) + 0.1 * span)
  n <- nrow(d)
  ypos <- rev(seq_len(n))
  mar <- graphics::par("mar")
  on.exit(graphics::par(mar = mar), add = TRUE)
  graphics::par(mar = c(4, max(mar[2L], 0.45 * max(nchar(d$function_name)) + 4),
                         2, 1))
  args <- utils::modifyList(
    list(xlab = "extinction time", ylab = "", yaxt = "n", xlim = xlim,
         ylim = c(0.5, n + 0.5), type = "n"),
    list(...)
  )
  do.call(graphics::plot, c(list(NA_real_, NA_real_), args))
  graphics::axis(2L, at = ypos, labels = paste0(d$function_name, "()"),
                 las = 1, cex.axis = 0.8, tick = FALSE)
  graphics::abline(h = ypos, col = "grey90")
  for (i in seq_len(n)) {
    lo <- d$lower[i]
    hi <- d$upper[i]
    if (!is.na(lo) || !is.na(hi)) {
      a <- if (is.na(lo)) NA_real_ else lo
      b <- if (is.na(hi)) NA_real_ else if (is.infinite(hi)) xlim[2L] else hi
      if (!is.na(a) && !is.na(b)) {
        graphics::segments(a, ypos[i], b, ypos[i], lwd = 3, col = "steelblue")
      }
    }
    if (!is.na(d$estimate[i])) graphics::points(d$estimate[i], ypos[i],
                                                pch = 19)
    if (is.na(d$estimate[i]) && is.finite(hi)) {
      graphics::points(hi, ypos[i], pch = 25, bg = "steelblue")
    }
  }
  invisible(x)
}

#' Coerce an extinction date estimate to a data frame
#'
#' @param x An [ede_estimate] object.
#' @param row.names,optional Not used; present for consistency with
#'   [base::as.data.frame()].
#' @param ... Not used.
#'
#' @return A one-row data frame with columns `method`, `estimate`, `lower`,
#'   `upper`, `interval_type` and `alpha`.
#'
#' @examples
#' as.data.frame(ole(sighting_data(dodo)))
#'
#' @export
as.data.frame.ede_estimate <- function(x, row.names = NULL, optional = FALSE,
                                       ...) {
  data.frame(
    method = x$method,
    estimate = x$estimate,
    lower = x$lower,
    upper = x$upper,
    interval_type = if (is.null(x$interval_type)) NA_character_
                    else x$interval_type,
    alpha = x$alpha,
    stringsAsFactors = FALSE,
    row.names = row.names
  )
}
