#' Thompson et al. (2013) probability of persistence with uncertain sightings
#'
#' Discrete-time model of Thompson, Lee, Stone, McCarthy & Burgman (2013) for
#' a record of one or more observation types (e.g. specimens, photographs,
#' sound recordings) seen or not seen in each time unit. The first type is
#' definite (no false sightings); the others may be false. For type `alpha`,
#' `p[alpha]` is the probability of not seeing that type in a time unit when
#' the species is extant and `q[alpha]` the probability of not seeing it when
#' the species is extinct (`q = 1` for the definite type). Time units are
#' independent, so the probability of the record given persistence to the
#' last time unit `T` is the product of the `p_i` (equation 5), and given
#' extinction in time unit `j` it is `p_1 ... p_(j-1) q_j ... q_T`
#' (equation 7), where `p_i` and `q_i` are products over types of `p` or
#' `1 - p` (`q` or `1 - q`) according to whether the type was seen
#' (equation 27). The probability of persistence is equation 9,
#'
#' `P(X_T | s) = prod(p_i) P(X_T) / (prod(p_i) P(X_T) + sum_j p_1...p_(j-1) q_j...q_T P_E)`,
#'
#' with `j` running over the time units after the last definite sighting
#' `T_N`, `P(X_T)` the prior probability of persistence and `P_E` the prior
#' probability of extinction in each of those time units.
#'
#' Each `p[alpha]` and `q[alpha]` is either a point value or a range
#' `c(a, b)` read as a uniform distribution, as in the paper. With point
#' values the result is the deterministic estimate `P_D` (equation 14 when
#' there is only the definite type). With ranges the result is averaged
#' (Section 2.2): `"annealed"` averages the numerator and the denominator of
#' equation 9 separately (a ratio of averages; equation 12), `"quenched"`
#' averages the ratio itself (equation 11). Parameters of different types,
#' and `p` and `q` of the same type, are independent. Each parameter is held
#' fixed across time units within a realisation.
#'
#' The record is binary: a time unit counts as a sighting of a type when its
#' count is positive. Time units are numbered from 1 at the earliest time of
#' any record (so, unlike [solow1993()], the first sighting is part of the
#' record), and `T = test_year - origin + 1`. Times must be whole numbers.
#'
#' Package constructs: the paper estimates the averages by simulation; here
#' the annealed averages are computed exactly, as products of
#' one-dimensional integrals ([stats::integrate()], on the log scale), and
#' the quenched average by [stats::integrate()] for a single random parameter
#' or by tensor-product Gauss-Legendre quadrature otherwise. The Bayes
#' factor (deterministic and annealed only) is the likelihood under
#' persistence divided by the mean likelihood over the extinction years,
#' which gives the paper's equation 13 in its two-year example; it is `NA`
#' for the quenched average, which is not a ratio of likelihoods.
#'
#' @param definite A [sighting_data] object for the definite (certain)
#'   observation type.
#' @param uncertain `NULL`, a [sighting_data] object, or a list of them, one
#'   per uncertain observation type.
#' @param test_year Last time unit `T` of the record. Must be later than the
#'   last definite sighting and not earlier than any other sighting.
#' @param p List (or numeric vector, for point values) with one element per
#'   observation type, definite first: the probability of not observing that
#'   type in a time unit given the species is extant. Each element is a point
#'   value in (0, 1) or a range `c(a, b)` with `0 <= a < b <= 1`.
#' @param q As `p`, for the uncertain types only: the probability of not
#'   observing that type given the species is extinct. Point values may be 1.
#' @param prior_extant Prior probability `P(X_T)` that the species is extant
#'   at `test_year`. Default `T_N / T` (equation 19).
#' @param prior_extinct Prior probability `P_E` of extinction in each time
#'   unit after the last definite sighting. Default `1 / T` (equation 19).
#' @param average `"annealed"` (default) or `"quenched"`; ignored when all
#'   parameters are point values.
#' @param nodes Number of Gauss-Legendre nodes per random parameter for the
#'   quenched average with more than one random parameter. Default chooses
#'   about 2e5 nodes in total.
#'
#' @return A list of class `ede_bayes` with components `bayes_factor` (in
#'   favour of persistence, or `NA`), `posterior` (probability of persistence
#'   from equation 9), `prior` (`prior_extant`), `test_year`, `method`, plus
#'   `prior_extinct` and `average` (`"deterministic"`, `"annealed"` or
#'   `"quenched"`).
#'
#' @references
#' Thompson, C. J., Lee, T. E., Stone, L., McCarthy, M. A. & Burgman, M. A.
#'   (2013). Inferring extinction risks from sighting records. Journal of
#'   Theoretical Biology, 338, 16-22.
#'
#' @export
thompson2013 <- function(definite, uncertain = NULL, test_year, p, q = NULL,
                         prior_extant = NULL, prior_extinct = NULL,
                         average = c("annealed", "quenched"), nodes = NULL) {
  average <- match.arg(average)
  if (missing(p)) p <- NULL
  if (!inherits(definite, "sighting_data")) {
    stop("`definite` must be a sighting_data object.", call. = FALSE)
  }
  if (inherits(uncertain, "sighting_data")) uncertain <- list(uncertain)
  if (is.null(uncertain)) uncertain <- list()
  if (!is.list(uncertain) ||
      !all(vapply(uncertain, inherits, logical(1), "sighting_data"))) {
    stop("`uncertain` must be NULL, a sighting_data object or a list of them.",
         call. = FALSE)
  }
  records <- c(list(definite), uncertain)
  n_type <- length(records)
  if (!any(definite$count > 0)) stop("no definite sightings with count > 0.", call. = FALSE)
  all_time <- unlist(lapply(records, function(r) r$time))
  if (any(all_time != round(all_time))) {
    stop("Thompson et al. (2013) needs whole-number time units.", call. = FALSE)
  }
  origin <- min(all_time)
  last_definite <- max(definite$time[definite$count > 0])
  last_any <- max(unlist(lapply(records, function(r) r$time[r$count > 0])))
  check_test_year(test_year, last_definite)
  if (test_year != round(test_year)) stop("`test_year` must be a whole number.", call. = FALSE)
  if (test_year < last_any) {
    stop("`test_year` cannot precede the last sighting.", call. = FALSE)
  }
  Tt <- test_year - origin + 1
  TN <- last_definite - origin + 1
  seen <- vapply(records, function(r) {
    s <- numeric(Tt)
    pos <- r$time[r$count > 0] - origin + 1
    pos <- pos[pos <= Tt]
    s[pos] <- 1
    s
  }, numeric(Tt))
  seen <- matrix(seen, nrow = Tt)

  p <- thompson2013_params(p, n_type, "p", allow_one = FALSE)
  q <- if (n_type > 1) thompson2013_params(q, n_type - 1L, "q", allow_one = TRUE) else list()
  if (is.null(prior_extant)) prior_extant <- TN / Tt
  if (is.null(prior_extinct)) prior_extinct <- 1 / Tt
  for (nm in c("prior_extant", "prior_extinct")) {
    v <- get(nm)
    if (length(v) != 1L || !is.numeric(v) || is.na(v) || v <= 0 || v >= 1) {
      stop("`", nm, "` must be in (0, 1).", call. = FALSE)
    }
  }
  # after the last definite sighting the definite type is never seen; q = 1
  after <- (TN + 1):Tt
  random <- any(lengths(c(p, q)) == 2L)
  if (!random) average <- "deterministic"

  # log E[prod_i x_i] for one parameter, x_i = x (not seen) or 1 - x (seen)
  log_moment <- function(par, n_seen, n_unseen) {
    lf <- function(x) uncertain_xlog(n_seen, 1 - x) + uncertain_xlog(n_unseen, x)
    if (length(par) == 1L) return(lf(par))
    thompson2013_log_mean(lf, par[1], par[2])
  }

  if (average != "quenched") {
    # numerator: prod over types of E[prod_{i <= T} p_i]
    log_num <- log(prior_extant) + sum(vapply(seq_len(n_type), function(a) {
      log_moment(p[[a]], sum(seen[, a]), Tt - sum(seen[, a]))
    }, numeric(1)))
    log_den <- vapply(after, function(j) {
      lp <- sum(vapply(seq_len(n_type), function(a) {
        s <- seen[seq_len(j - 1L), a]
        log_moment(p[[a]], sum(s), length(s) - sum(s))
      }, numeric(1)))
      lq <- sum(vapply(seq_len(n_type - 1L), function(a) {
        s <- seen[j:Tt, a + 1L]
        log_moment(q[[a]], sum(s), length(s) - sum(s))
      }, numeric(1)))
      log(prior_extinct) + lp + lq
    }, numeric(1))
    log_l_ext <- log_num - log(prior_extant)
    log_l_extinct <- thompson2013_logsumexp(log_den - log(prior_extinct)) -
      log(length(after))
    bf <- exp(log_l_ext - log_l_extinct)
    posterior <- 1 / (1 + exp(thompson2013_logsumexp(log_den) - log_num))
  } else {
    pars <- c(p, q)
    is_rand <- lengths(pars) == 2L
    post_at <- function(vals) {
      # vals: matrix, one row per realisation, one column per parameter
      pp <- vals[, seq_len(n_type), drop = FALSE]
      qq <- vals[, n_type + seq_len(n_type - 1L), drop = FALSE]
      thompson2013_posterior(seen, pp, qq, TN, prior_extant, prior_extinct)
    }
    fixed <- vapply(pars, function(x) x[1], numeric(1))
    if (sum(is_rand) == 1L) {
      r <- pars[[which(is_rand)]]
      f <- function(x) {
        vals <- matrix(fixed, length(x), length(fixed), byrow = TRUE)
        vals[, which(is_rand)] <- x
        post_at(vals)
      }
      posterior <- stats::integrate(f, r[1], r[2], rel.tol = 1e-10,
                                    subdivisions = 1000L)$value / (r[2] - r[1])
    } else {
      d <- sum(is_rand)
      if (is.null(nodes)) nodes <- max(4L, floor(2e5^(1 / d)))
      gl <- thompson2013_gauss_legendre(nodes)
      axes <- lapply(pars[is_rand], function(r) (r[1] + r[2]) / 2 + (r[2] - r[1]) / 2 * gl$x)
      grid <- as.matrix(expand.grid(axes))
      w <- Reduce(`*`, expand.grid(rep(list(gl$w / 2), d)))
      vals <- matrix(fixed, nrow(grid), length(fixed), byrow = TRUE)
      vals[, is_rand] <- grid
      posterior <- sum(w * post_at(vals))
    }
    bf <- NA_real_
  }
  structure(
    list(bayes_factor = bf, posterior = posterior, prior = prior_extant,
         test_year = test_year,
         method = sprintf("Thompson et al. (2013) %s", average),
         prior_extinct = prior_extinct, average = average),
    class = "ede_bayes"
  )
}

#' Validate a p or q specification into a list of points / ranges
#' @keywords internal
#' @noRd
thompson2013_params <- function(x, n, name, allow_one) {
  if (is.null(x)) stop("`", name, "` must be supplied.", call. = FALSE)
  if (is.numeric(x)) x <- as.list(x)
  if (!is.list(x) || length(x) != n) {
    stop("`", name, "` needs one element per ",
         if (name == "p") "observation type" else "uncertain observation type",
         " (", n, ").", call. = FALSE)
  }
  lapply(x, function(v) {
    ok <- is.numeric(v) && !anyNA(v) && length(v) %in% 1:2
    if (ok && length(v) == 1L) {
      ok <- v > 0 && (v < 1 || (allow_one && v == 1))
    } else if (ok) {
      ok <- v[1] >= 0 && v[2] <= 1 && v[1] < v[2]
    }
    if (!ok) {
      stop("each element of `", name, "` must be a probability in (0, 1",
           if (allow_one) "]" else ")", " or a range c(a, b) with ",
           "0 <= a < b <= 1.", call. = FALSE)
    }
    v
  })
}

#' Posterior of persistence (eq 9) for each row of parameter values
#' @keywords internal
#' @noRd
thompson2013_posterior <- function(seen, pp, qq, TN, prior_extant, prior_extinct) {
  Tt <- nrow(seen)
  G <- nrow(pp)
  logp <- matrix(0, G, Tt)
  for (a in seq_len(ncol(seen))) {
    s <- seen[, a] == 1
    lp <- matrix(log(pp[, a]), G, Tt)
    lp[, s] <- log1p(-pp[, a])
    logp <- logp + lp
  }
  after <- (TN + 1):Tt
  logq <- matrix(0, G, length(after))
  for (a in seq_len(ncol(qq))) {
    s <- seen[after, a + 1L] == 1
    lq <- matrix(log(qq[, a]), G, length(after))
    lq[, s] <- log1p(-qq[, a])
    logq <- logq + lq
  }
  cum_p <- matrix(0, G, Tt + 1L)                    # cum_p[, j] = sum_{i < j}
  for (j in seq_len(Tt)) cum_p[, j + 1L] <- cum_p[, j] + logp[, j]
  n_after <- length(after)
  tail_q <- matrix(0, G, n_after)                   # sum_{i >= j}, i, j > TN
  tail_q[, n_after] <- logq[, n_after]
  for (j in rev(seq_len(n_after - 1L))) tail_q[, j] <- tail_q[, j + 1L] + logq[, j]
  log_num <- log(prior_extant) + cum_p[, Tt + 1L]
  log_den <- log(prior_extinct) + cum_p[, after, drop = FALSE] + tail_q
  m <- apply(log_den, 1, max)
  lse <- m + log(rowSums(exp(log_den - m)))
  lse[!is.finite(m)] <- -Inf
  1 / (1 + exp(lse - log_num))
}

#' log of the mean of exp(lf(x)) for x uniform on (a, b)
#' @keywords internal
#' @noRd
thompson2013_log_mean <- function(lf, a, b) {
  grid <- seq(a, b, length.out = 401L)
  inner <- grid[-c(1L, length(grid))]
  vals <- lf(inner)
  m <- max(vals[is.finite(vals)])
  opt <- stats::optimize(lf, c(a, b), maximum = TRUE)
  if (is.finite(opt$objective) && opt$objective > m) m <- opt$objective
  mid <- if (is.finite(opt$objective) && opt$objective >= max(vals)) opt$maximum else inner[which.max(vals)]
  g <- function(x) exp(lf(x) - m)
  v <- stats::integrate(g, a, mid, rel.tol = 1e-10, subdivisions = 1000L)$value +
    stats::integrate(g, mid, b, rel.tol = 1e-10, subdivisions = 1000L)$value
  m + log(v / (b - a))
}

#' @keywords internal
#' @noRd
thompson2013_logsumexp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

#' Gauss-Legendre nodes and weights on (-1, 1) (Golub-Welsch)
#' @keywords internal
#' @noRd
thompson2013_gauss_legendre <- function(n) {
  i <- seq_len(n - 1L)
  b <- i / sqrt(4 * i^2 - 1)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1L)] <- b
  J[cbind(i + 1L, i)] <- b
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = 2 * e$vectors[1, ]^2)
}
