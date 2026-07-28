#' Burgman, Grimson & Ferson (1995) combinatorial persistence test
#'
#' Computes the probability that, if sighting events were distributed
#' uniformly at random over the candidate observation window, the largest
#' observed gap between sightings would not exceed the gap actually seen in
#' the data. Uses an inclusion-exclusion (Stirling-number) argument on the
#' occupancy of time bins by sighting events.
#'
#' @inheritParams solow1993
#'
#' @return An [ede_estimate] object, or (if `data_out = TRUE`) a data frame
#'   with columns `time` and `chance`.
#'
#' @references
#' Burgman, M. A., Grimson, R. C., & Ferson, S. (1995). Inferring threat
#'   from scientific collections. Conservation Biology, 9(4), 923-928.
#'
#' @export
burgman1995 <- function(sd, alpha = 0.05, test_year, data_out = FALSE) {
  stopifnot(inherits(sd, "sighting_data"))
  if (alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0, 1).", call. = FALSE)
  if (missing(test_year) || !is.numeric(test_year)) {
    stop("`test_year` must be supplied as a number.", call. = FALSE)
  }

  full <- expand_record(sd, test_year)
  last_sight <- max(sd$time[sd$count > 0])
  candidates <- full$time[full$time > last_sight]
  if (length(candidates) == 0L) {
    stop("`test_year` must be later than the last sighting.", call. = FALSE)
  }

  chance <- vapply(candidates, function(t) {
    burgman_chance(full[full$time <= t, , drop = FALSE])
  }, numeric(1))

  if (data_out) return(data.frame(time = candidates, chance = chance))

  decreasing <- c(FALSE, diff(chance) < 0)
  below <- candidates[decreasing & chance <= alpha]
  if (length(below) == 0L) {
    warning(
      "chance of persistence never falls to alpha (on a decreasing run) ",
      "before `test_year`; returning NA.", call. = FALSE
    )
    return(new_ede_estimate(NA_real_, method = "Burgman, Grimson & Ferson (1995)", alpha = alpha))
  }
  new_ede_estimate(below[1], method = "Burgman, Grimson & Ferson (1995)", alpha = alpha)
}

#' @keywords internal
#' @noRd
stirling2_matrix <- function(n_max) {
  S <- matrix(0, nrow = n_max, ncol = n_max)
  S[1, 1] <- 1
  if (n_max > 1) {
    for (n in 2:n_max) {
      S[n, 1] <- 1
      for (j in 2:n) {
        S[n, j] <- j * S[n - 1, j] + S[n - 1, j - 1]
      }
    }
  }
  S
}

#' @keywords internal
#' @noRd
burgman_equation4_pvalue <- function(ct, n, r) {
  if (n <= 0 || r <= 0 || r >= ct) return(1.0)

  s_mat <- stirling2_matrix(n)
  total_sum <- 0

  for (j in seq_len(min(n, ct))) {
    s_nj <- s_mat[n, j]
    if (s_nj == 0) next

    k_max <- min(j + 1, floor((ct - j) / r))
    if (k_max < 1) next

    for (k in seq_len(k_max)) {
      sign_k <- if (k %% 2 == 1) 1 else -1

      a <- ct - r * k
      log_fallfac <- lgamma(a + 1) - lgamma(a - j + 1)
      log_comb <- lchoose(j + 1, k)

      log_term <- log_comb + log_fallfac + log(s_nj)
      term <- sign_k * exp(log_term)

      total_sum <- total_sum + term
    }
  }

  p <- total_sum * (ct^(-n))
  p <- max(0, min(1, p))
  p
}

#' @keywords internal
#' @noRd
burgman_chance <- function(d) {
  time_vec <- d$time
  count_vec <- d$count

  n <- sum(count_vec)
  ct <- length(time_vec)

  if (n == 0 || ct == 0) return(1.0)

  nz_indices <- which(count_vec > 0)
  if (length(nz_indices) == 0) return(1.0)

  gap_start <- nz_indices[1] - 1
  gaps_mid  <- if (length(nz_indices) > 1) diff(nz_indices) - 1 else 0
  gap_end   <- ct - nz_indices[length(nz_indices)]

  r <- max(c(gap_start, gaps_mid, gap_end))

  burgman_equation4_pvalue(ct, n, r)
}
