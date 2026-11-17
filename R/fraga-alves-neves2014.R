#' Fraga Alves & Neves (2014) end-point estimator for the Gumbel domain
#'
#' Estimates the extinction time as the finite right end-point of the sighting
#' time distribution when that distribution lies in the Gumbel domain of
#' attraction: a density that vanishes at the end-point faster than any power
#' of the distance to it (for example an exponentially declining sighting rate
#' that is cut off at the extinction time). It complements [ole()] and
#' [hall_wang1999()], which assume the Weibull domain, a density that behaves
#' like a power of the distance to the end-point; the Weibull models become
#' increasingly Gumbel-like as their extremal exponent grows.
#'
#' With the sighting times ordered `t[1] > t[2] > ...` (latest first) the
#' estimate is
#' `t[1] + sum_{i = 0}^{k - 1} a[i] * (t[k + 1] - t[k + i + 1])`, with weights
#' `a[i] = (log(k + i + 1) - log(k + i)) / log(2)` that are positive and sum
#' to one. It is therefore never earlier than the last sighting, which the
#' usual estimators for the Weibull domain do not guarantee, and it uses the
#' `2k` latest sightings.
#'
#' The limiting distribution of the estimator depends on unknown auxiliary
#' functions of the sampling distribution, so Fraga Alves & Neves (2014)
#' provide no interval, and neither does this function: `lower` and `upper`
#' are `NA`.
#'
#' Each distinct sighting time counts once by default, as in [robson1964()];
#' set `events = TRUE` to count every sighting event, as [ole()] does.
#'
#' @param sd A [sighting_data] object.
#' @param k Number of top order statistics defining the high threshold, a
#'   positive integer with `2k` at most the number of sightings, or `NULL`
#'   for the largest admissible value (`floor(n / 2)`).
#' @param events If `FALSE` (default), count each distinct sighting time once;
#'   if `TRUE`, count every sighting event.
#'
#' @return An [ede_estimate] object with `estimate` and no interval, plus the
#'   component `k`.
#'
#' @references
#' Fraga Alves, I., & Neves, C. (2014). Estimation of the finite right
#'   endpoint in the Gumbel domain. Statistica Sinica, 24(4), 1811-1835.
#'   doi:10.5705/ss.2013.183
#'
#' @seealso [hall_wang1999()] and [ole()] for the Weibull domain.
#'
#' @examples
#' sd <- sighting_data(data.frame(year = c(1900, 1904, 1909, 1911, 1913, 1918,
#'                                          1920, 1925, 1928, 1930, 1934),
#'                                sightings = 1))
#' fraga_alves_neves2014(sd)
#' fraga_alves_neves2014(sd, k = 3)
#'
#' @export
fraga_alves_neves2014 <- function(sd, k = NULL, events = FALSE) {
  stopifnot(inherits(sd, "sighting_data"))
  t <- evt_times(sd, events)
  n <- length(t)
  if (n < 2L) {
    stop("Fraga Alves & Neves (2014) needs at least 2 sighting times; the ",
         "record has ", n, ".", call. = FALSE)
  }
  k <- evt_k(k, n, kmin = 1L, default = n %/% 2L,
             label = "Fraga Alves & Neves (2014)", kmax = n %/% 2L)
  i <- 0:(k - 1L)
  a <- (log(k + i + 1L) - log(k + i)) / log(2)
  estimate <- t[1L] + sum(a * (t[k + 1L] - t[k + i + 1L]))
  out <- new_ede_estimate(estimate, method = "Fraga Alves & Neves (2014)",
                          alpha = NA_real_, interval_type = "two-sided")
  out$k <- k
  out
}
