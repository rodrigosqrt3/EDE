#' Burgman, Grimson & Ferson (1995) tests for frequency sighting data
#'
#' Exact tests for sighting records divided into `CT` equal periods, where `N`
#' is the total number of sightings (all records counted, not only occupied
#' periods).
#'
#' * `equation = 2` (the default) tests extinction. The p-value
#'   `(Ce / CT)^N`, where `Ce` is the period of the last sighting, is the
#'   probability that the species is not extinct; it is the empty-cells test
#'   restricted to the run of empty periods at the end of the record.
#' * `equation = 3` is the runs test of Bradley (1968) on the binary
#'   (occupied / empty) sequence of periods, so counts are ignored. With `n1`
#'   occupied periods out of `CT` and `r` the longest run of empty periods
#'   anywhere in the record, the probability of an empty run at least `r`
#'   long is
#'   `sum_{k=1}^{[n0/r]} (-1)^(k+1) C(n1+1, k) C(CT - r k, n1) / C(CT, n1)`,
#'   with `n0 = CT - n1` empty periods (inclusion-exclusion over the `n1 + 1`
#'   gaps that hold the empty periods). The equation is illegible in the
#'   source, so this form was reconstructed from its definitions and checked
#'   against the paper's numbers: it gives the Fig. 1b value 0.02885 and the
#'   Table 2 values for the black-footed ferret (0.003) and both stick-nest
#'   rats (0.139, 0.350). It does not give the Table 2 value for the monk seal
#'   (0.244 here against 0.258 in the paper, which equals the paper's
#'   equation 4 value; see Details). Like equation 4, it addresses a change in
#'   status rather than extinction.
#' * `equation = 4` tests for a change in conservation status (Grimson et al.,
#'   1992). It computes, by inclusion-exclusion over the occupancy of periods,
#'   the probability of a run of empty periods at least as long as the longest
#'   observed run anywhere in the record. Burgman et al. (1995, p. 927) note
#'   that it addresses a different question from extinction, and that a
#'   significant result may be discounted if the longest run is not near the
#'   end of the record.
#'
#' As in the Caribbean monk seal example of the source paper (Fig. 1a and
#' Tables 1-2), the earliest supplied time defines the observation origin and
#' is not one of the `CT` periods; sightings recorded in that period are not
#' counted. To include an earlier start date, prepend a zero-count row. The
#' paper's black-footed ferret values (observation from January 1972)
#' correspond to a zero-count row for 1971, so that 1972 is a period. Burgman
#' et al. (2000, Table 2.1) used the same convention for Western Australian
#' *Acacia*: the year of the first collection is the first period and 1993
#' the last; equations 2 and 4 reproduce their published values.
#'
#' For the monk seal, the origin convention above reproduces the paper's
#' equation 2 (0.053) and equation 4 (0.258) values (`CT = 77`, `n1 = N = 4`,
#' `r = 40`), but equation 3 gives 0.244. No choice of `CT` and `r` near
#' these values gives 0.258 exactly (`CT = 78` with `n1 = 4`, i.e. counting
#' 1915 as an empty period, gives 0.2587). Since all seal counts are 1, the
#' paper's equal entries for equations 3 and 4 suggest that the equation 4
#' value was reported for both.
#'
#' The paper defines only the tests. Since the periods are discrete, the
#' returned upper bound is the first period after the last sighting at which
#' the null hypothesis is rejected. Under `equation = 3` or `4` the p-value
#' curve need not be monotone, so this bound should be read with the caveat
#' above.
#'
#' @inheritParams solow1993
#' @param equation Which equation of the source paper to use: `2` for
#'   extinction (the default), `3` (binary runs test) or `4` (frequency
#'   runs test) for a change in conservation status.
#'
#' @return An [ede_estimate] object with `estimate = NA` and the one-sided
#'   upper bound in `upper`, or (if `data_out = TRUE`) a data frame with
#'   columns `time` and `chance` (p-values).
#'
#' @references
#' Burgman, M. A., Grimson, R. C., & Ferson, S. (1995). Inferring threat
#'   from scientific collections. Conservation Biology, 9(4), 923-928.
#'
#' Burgman, M. A., Maslin, B. R., Andrewartha, D., Keatley, M. R., Boek, C.,
#'   & McCarthy, M. (2000). Inferring threat from scientific collections:
#'   power tests and an application to Western Australian Acacia species. In
#'   S. Ferson & M. Burgman (Eds.), Quantitative Methods for Conservation
#'   Biology (pp. 7-26). Springer, New York.
#'
#' Bradley, J. V. (1968). Distribution-free statistical tests. Prentice-Hall,
#'   Englewood Cliffs, New Jersey.
#'
#' @export
burgman1995 <- function(sd, alpha = 0.05, test_year, data_out = FALSE,
                        equation = 2) {
  stopifnot(inherits(sd, "sighting_data"))
  if (alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0, 1).", call. = FALSE)
  if (length(equation) != 1L || !equation %in% c(2, 3, 4)) {
    stop("`equation` must be 2, 3 or 4.", call. = FALSE)
  }
  check_test_year(test_year)
  if (any(sd$time != floor(sd$time))) {
    stop("Burgman et al. (1995) needs whole-number times (one row per period).", call. = FALSE)
  }

  full <- expand_record(sd, test_year)
  last_sight <- max(sd$time[sd$count > 0])
  check_test_year(test_year, last_sight)
  candidates <- full$time[full$time > last_sight]

  chance_fun <- switch(as.character(equation),
                       "2" = burgman_equation2_chance,
                       "3" = burgman_equation3_chance,
                       "4" = burgman_chance)
  chance <- vapply(candidates, function(t) {
    chance_fun(full[full$time <= t, , drop = FALSE])
  }, numeric(1))

  if (data_out) return(data.frame(time = candidates, chance = chance))

  method <- sprintf("Burgman, Grimson & Ferson (1995), equation %d", equation)
  below <- candidates[chance <= alpha]
  bound <- if (length(below) > 0L) below[1] else NA_real_
  test_bound_result(bound, test_year, method, alpha)
}

#' Periods after the observation origin, as used by Burgman et al. (1995)
#' @keywords internal
#' @noRd
burgman_periods <- function(d) {
  d[d$time > min(d$time), , drop = FALSE]
}

#' @keywords internal
#' @noRd
burgman_equation2_chance <- function(d) {
  d <- burgman_periods(d)
  n <- sum(d$count)
  if (n == 0) return(1.0)
  ce <- max(which(d$count > 0))
  (ce / nrow(d))^n
}

#' Longest run of empty periods anywhere in a count vector
#' @keywords internal
#' @noRd
burgman_longest_empty_run <- function(count_vec) {
  runs <- rle(count_vec == 0)
  max(c(0, runs$lengths[runs$values]))
}

#' Equation 3: Bradley's runs test on the binary sequence of periods
#' @keywords internal
#' @noRd
burgman_equation3_pvalue <- function(ct, n1, r) {
  if (n1 <= 0 || r <= 0) return(1.0)
  if (r > ct - n1) return(0)
  k <- seq_len(min(n1 + 1, floor((ct - n1) / r)))
  log_terms <- lchoose(n1 + 1, k) + lchoose(ct - r * k, n1) - lchoose(ct, n1)
  signs <- ifelse(k %% 2 == 1, 1, -1)
  m <- max(log_terms)
  scaled <- signs * exp(log_terms - m)
  total <- sum(scaled)
  if (abs(total) < 1e-10 * sum(abs(scaled))) {
    warning("equation 3 lost precision through cancellation; the p-value is unreliable.",
            call. = FALSE)
  }
  p <- sign(total) * exp(m + log(abs(total)))
  max(0, min(1, p))
}

#' @keywords internal
#' @noRd
burgman_equation3_chance <- function(d) {
  d <- burgman_periods(d)
  if (nrow(d) == 0L || !any(d$count > 0)) return(1.0)
  burgman_equation3_pvalue(nrow(d), sum(d$count > 0),
                           burgman_longest_empty_run(d$count))
}

#' Second-kind Stirling numbers on the log scale
#' @keywords internal
#' @noRd
log_stirling2_matrix <- function(n_max) {
  S <- matrix(-Inf, nrow = n_max, ncol = n_max)
  S[, 1] <- 0
  if (n_max > 1) {
    for (n in 2:n_max) {
      for (j in 2:n) {
        a <- log(j) + S[n - 1, j]
        b <- S[n - 1, j - 1]
        m <- max(a, b)
        S[n, j] <- m + log(exp(a - m) + exp(b - m))
      }
    }
  }
  S
}

#' @keywords internal
#' @noRd
burgman_equation4_pvalue <- function(ct, n, r) {
  if (n <= 0 || r <= 0 || r >= ct) return(1.0)

  log_s <- log_stirling2_matrix(n)[n, ]
  log_terms <- numeric(0)
  signs <- numeric(0)

  for (j in seq_len(min(n, ct))) {
    k_max <- min(j + 1, floor((ct - j) / r))
    if (k_max < 1) next
    k <- seq_len(k_max)
    a <- ct - r * k
    log_terms <- c(log_terms, lchoose(j + 1, k) + lgamma(a + 1) -
                     lgamma(a - j + 1) + log_s[j] - n * log(ct))
    signs <- c(signs, ifelse(k %% 2 == 1, 1, -1))
  }
  if (length(log_terms) == 0L) return(0)

  # scale by the largest term so neither the Stirling numbers nor ct^n
  # overflow; cancellation between alternating terms is checked below
  m <- max(log_terms)
  scaled <- signs * exp(log_terms - m)
  total <- sum(scaled)
  if (abs(total) < 1e-10 * sum(abs(scaled))) {
    warning("equation 4 lost precision through cancellation; the p-value is unreliable.",
            call. = FALSE)
  }
  p <- sign(total) * exp(m + log(abs(total)))
  max(0, min(1, p))
}

#' @keywords internal
#' @noRd
burgman_chance <- function(d) {
  d <- burgman_periods(d)
  count_vec <- d$count

  n <- sum(count_vec)
  ct <- length(count_vec)

  if (n == 0 || ct == 0) return(1.0)

  nz_indices <- which(count_vec > 0)

  gap_start <- nz_indices[1] - 1
  gaps_mid  <- if (length(nz_indices) > 1) diff(nz_indices) - 1 else 0
  gap_end   <- ct - nz_indices[length(nz_indices)]

  r <- max(c(gap_start, gaps_mid, gap_end))

  burgman_equation4_pvalue(ct, n, r)
}
