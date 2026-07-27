#' McInerny, Roberts, Davy & Cribb (2006) sighting-rate persistence test
#'
#' A modification of the Solow (1993) persistence test that conditions on
#' the sighting rate observed up to the last sighting, `n / tn`, instead of
#' on the length of the whole observation window. This makes the test
#' comparable across records with very different total observation
#' periods: species discovered recently and species discovered long ago,
#' but sighted at the same rate, are inferred extinct after the same
#' length of silence.
#'
#' Following the source paper, the first sighting is used to anchor the
#' time origin (`t = 0`) rather than counted as one of the `n` sighting
#' events being tested, so the count entering the formula is `n - 1`,
#' where `n` is the number of distinct sighting times with a positive
#' count.
#'
#' @param sd A [sighting_data] object.
#' @param alpha Significance level, in (0, 1). Persistence is rejected for
#'   the first candidate year at which the chance of persistence falls to
#'   or below `alpha`.
#' @param test_year Latest year to test. Must be later than the last
#'   sighting.
#' @param data_out If `TRUE`, return the full chance-of-persistence curve
#'   instead of the single first-rejection year.
#'
#' @return An [ede_estimate] object, or (if `data_out = TRUE`) a data frame
#'   with columns `time` and `chance`.
#'
#' @references
#' McInerny, G. J., Roberts, D. L., Davy, A. J., & Cribb, P. J. (2006).
#'   Significance of sighting rate in inferring extinction and threat.
#'   Conservation Biology, 20(2), 562-567.
#'
#' @export
mcinerny2006 <- function(sd, alpha = 0.05, test_year, data_out = FALSE) {
  stopifnot(inherits(sd, "sighting_data"))
  if (alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0, 1).", call. = FALSE)
  if (missing(test_year) || !is.numeric(test_year)) {
    stop("`test_year` must be supplied as a number.", call. = FALSE)
  }

  times <- sort(unique(sd$time[sd$count > 0]))
  n <- length(times)
  if (n < 2) {
    stop(
      "McInerny et al. (2006) needs at least 2 distinct sighting times ",
      "with count > 0.", call. = FALSE
    )
  }

  last_sight <- times[n]
  if (test_year <= last_sight) {
    stop("`test_year` must be later than the last sighting.", call. = FALSE)
  }

  # first sighting anchors t = 0 and is not itself one of the n tested
  # sighting events (McInerny et al. 2006, Methods)
  t1 <- times[1]
  tn <- last_sight - t1
  n_used <- n - 1
  rate <- n_used / tn
  if (rate >= 1) {
    warning(
      "estimated sighting rate n / tn >= 1 (a sighting in essentially ",
      "every time unit); the sighting-rate equation is degenerate for ",
      "this record and will reject persistence immediately.", call. = FALSE
    )
  }

  candidates <- seq(last_sight + 1, test_year)
  chance <- (1 - rate)^((candidates - t1) - tn)

  if (data_out) return(data.frame(time = candidates, chance = chance))

  below <- candidates[chance <= alpha]
  if (length(below) == 0L) {
    warning(
      "chance of persistence never falls to alpha before `test_year`; ",
      "returning NA.", call. = FALSE
    )
    return(new_ede_estimate(NA_real_, method = "McInerny, Roberts, Davy & Cribb (2006)", alpha = alpha))
  }
  new_ede_estimate(below[1], method = "McInerny, Roberts, Davy & Cribb (2006)", alpha = alpha)
}
