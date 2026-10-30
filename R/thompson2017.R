#' Thompson et al. (2017) iterative probability that a taxon is extant
#'
#' Year-by-year model of the probability `P(X[t])` that a taxon is extant,
#' combining records of varying reliability with unsuccessful dedicated and
#' passive surveys (Thompson et al., 2017; the model behind the IUCN Red List
#' extinction-probability tool). Starting from an initial probability
#' `p_initial` at the beginning of `start_year`, each year updates the
#' probability by equation 5 of the paper:
#'
#' * record year: `P(X[t+1]) = p(ci) + (1 - p(ci)) * P(X[t])` (equation 1),
#'   where `p(ci)` is the probability that the record is a correct
#'   identification of the extant taxon;
#' * unsuccessful survey year: `P(X[t+1]) = (1 - eps * p(r) * p(i)) * P(X[t])`
#'   (equations 3-4), where `eps` is the proportion of the habitat surveyed,
#'   `p(r)` the probability that the taxon would have been recorded and
#'   `p(i)` the probability that it would have been reliably identified;
#' * gap year (neither a record nor a dedicated survey): the survey rule with
#'   the passive-survey inputs `eps'`, `p'(r)` and `p'(i)` (section 4.2).
#'
#' The model is not Bayesian (section 3 of the paper): `P(X[t])` is an
#' absolute probability updated by the rules above. Successful surveys are
#' supplied as records.
#'
#' Every input is given as a lower and upper bound, as in Table 1 of the
#' paper, and uncertainty is propagated in the two ways of section 6:
#'
#' * `p_extant` iterates equation 5 with the midpoint of every bound (the
#'   solid line of Figs. 1-3);
#' * `lower` and `upper` are the results of interval arithmetic. Equation 5
#'   increases with `P(X[t])` and `p(ci)` and decreases with `eps`, `p(r)`
#'   and `p(i)`, so the exact bounds are obtained by iterating it with the
#'   corresponding extreme inputs (the light shaded region of Figs. 1-3);
#' * if `n_sim > 0`, each bound pair is also treated as the 5th and 95th
#'   percentiles of an independent normal distribution, inputs are drawn
#'   independently for every year (including each gap year) and every
#'   replicate, and the 5th and 95th percentiles of the simulated
#'   `P(X[t])` are returned as `mc_lower` and `mc_upper` (the dark shaded
#'   region). The paper does not say how draws outside `[0, 1]` are handled;
#'   here they are truncated to `[0, 1]` (a package construct).
#'
#' A point value can be supplied instead of a range by using a single column
#' (e.g. `pci` instead of `pci_lower`/`pci_upper`).
#'
#' @param records Data frame of record years with columns `year` and
#'   `pci_lower`, `pci_upper` (or a single column `pci`), the probability of
#'   correct identification. `NULL` if there are no records.
#' @param surveys Data frame of unsuccessful dedicated survey years with
#'   columns `year`, `eps_lower`, `eps_upper`, `pi_lower`, `pi_upper`,
#'   `pr_lower`, `pr_upper` (or single columns `eps`, `pi`, `pr`). `NULL` if
#'   there are no dedicated surveys.
#' @param passive Passive-survey inputs applied to every gap year: a named
#'   numeric vector or list with elements `eps_lower`, `eps_upper`,
#'   `pi_lower`, `pi_upper`, `pr_lower`, `pr_upper` (or `eps`, `pi`, `pr`).
#'   Use `c(eps = 0, pi = 0, pr = 0)` for no passive surveillance.
#' @param p_initial Probability that the taxon is extant at the beginning of
#'   `start_year`: one value or a `c(lower, upper)` range.
#' @param start_year First year of the iteration; `p_initial` applies at its
#'   beginning. Defaults to the earliest year in `records` and `surveys`.
#' @param end_year Last year of the iteration. Defaults to the latest year in
#'   `records` and `surveys`; years after the last record or survey are gap
#'   years. Records or surveys outside `start_year:end_year` are ignored.
#'   Both years are required when there are no records or surveys (e.g. the
#'   passive-only projection of section 8.2 and equation 14 of the paper).
#' @param n_sim Number of Monte Carlo replicates for the percentile band; `0`
#'   skips the simulation.
#' @param data_out If `TRUE`, return the year-by-year trajectory instead of
#'   the summary object.
#'
#' @return A list of class `ede_extant` with components `p_extant`
#'   (midpoint estimate of the probability that the taxon is extant at the
#'   end of `end_year`), `lower`, `upper` (interval-arithmetic bounds),
#'   `mc_lower`, `mc_upper` (Monte Carlo 5th and 95th percentiles, or `NA`
#'   when `n_sim = 0`), `p_extinct` (`1 - p_extant`), `end_year` and
#'   `method`; or (if `data_out = TRUE`) a data frame with one row per year
#'   and columns `time`, `type` (`"record"`, `"survey"` or `"passive"`),
#'   `p_extant`, `lower`, `upper`, `mc_lower` and `mc_upper`, where each row
#'   gives the probability that the taxon is extant at the end of year
#'   `time`, after that year's data (as plotted in Figs. 1-3).
#'
#' @references
#' Thompson, C. J., Koshkina, V., Burgman, M. A., Butchart, S. H. M., &
#'   Stone, L. (2017). Inferring extinctions II: a practical, iterative model
#'   based on records and surveys. Biological Conservation, 214, 328-335.
#'
#' @examples
#' # Alaotra Grebe, Table 1 of Thompson et al. (2017)
#' rec <- data.frame(
#'   year = c(1929, 1960, 1963, 1969, 1970, 1971, 1972, 1982, 1985, 1986, 1988),
#'   pci_lower = c(0.99, 0.95, 0.75, 0.75, 0.10, 0.10, 0.60, 0.60, 0.20, 0.20, 0.20),
#'   pci_upper = c(1.00, 0.99, 0.94, 0.94, 0.40, 0.40, 0.80, 0.80, 0.80, 0.70, 0.50)
#' )
#' sv_years <- c(1989, 1990, 1993, 1994, 1997, 1998, 1999, 2000, 2004, 2009)
#' sv <- data.frame(
#'   year = sv_years,
#'   eps_lower = ifelse(sv_years == 2000, 0.70, 0.80),
#'   eps_upper = ifelse(sv_years == 2000, 0.90, 0.95),
#'   pi_lower = 0.90, pi_upper = 0.95, pr_lower = 0.70, pr_upper = 0.90
#' )
#' pas <- c(eps_lower = 0, eps_upper = 0.05, pi_lower = 0.10, pi_upper = 0.65,
#'          pr_lower = 0.40, pr_upper = 0.60)
#' thompson2017(rec, sv, pas, p_initial = c(0.95, 0.99), end_year = 2017,
#'              n_sim = 0)
#'
#' @export
thompson2017 <- function(records, surveys, passive, p_initial,
                         start_year = NULL, end_year = NULL,
                         n_sim = 10000, data_out = FALSE) {
  inp <- thompson_inputs(records, surveys, passive, p_initial,
                         start_year, end_year)
  if (length(n_sim) != 1L || !is.numeric(n_sim) || is.na(n_sim) ||
      n_sim < 0 || n_sim != round(n_sim)) {
    stop("`n_sim` must be a non-negative whole number.", call. = FALSE)
  }

  # lo/hi: bounds of each input; the extant probability is increasing in
  # P(X) and p(ci) and decreasing in the survey inputs (equation 5)
  run <- function(p0, pci, eps, pi, pr) {
    p <- numeric(inp$n)
    cur <- p0
    for (j in seq_len(inp$n)) {
      cur <- if (inp$type[j] == "record") {
        pci[j] + (1 - pci[j]) * cur                    # equation 1
      } else {
        (1 - eps[j] * pr[j] * pi[j]) * cur             # equations 3-4
      }
      p[j] <- cur
    }
    p
  }
  mid <- function(x) (x$lower + x$upper) / 2
  p_mid <- run(mid(inp$p0), mid(inp$pci), mid(inp$eps), mid(inp$pi), mid(inp$pr))
  p_low <- run(inp$p0$lower, inp$pci$lower, inp$eps$upper, inp$pi$upper,
               inp$pr$upper)
  p_upp <- run(inp$p0$upper, inp$pci$upper, inp$eps$lower, inp$pi$lower,
               inp$pr$lower)

  mc_lower <- mc_upper <- rep(NA_real_, inp$n)
  if (n_sim > 0) {
    z95 <- stats::qnorm(0.95)
    draw <- function(x, size) {
      m <- (x$lower + x$upper) / 2
      s <- (x$upper - x$lower) / (2 * z95)
      v <- stats::rnorm(size * length(m), mean = rep(m, each = size),
                        sd = rep(s, each = size))
      matrix(pmin(pmax(v, 0), 1), nrow = size)
    }
    cur <- as.vector(draw(inp$p0, n_sim))
    pci <- draw(inp$pci, n_sim)
    eps <- draw(inp$eps, n_sim)
    pi <- draw(inp$pi, n_sim)
    pr <- draw(inp$pr, n_sim)
    for (j in seq_len(inp$n)) {
      cur <- if (inp$type[j] == "record") {
        pci[, j] + (1 - pci[, j]) * cur
      } else {
        (1 - eps[, j] * pr[, j] * pi[, j]) * cur
      }
      q <- stats::quantile(cur, c(0.05, 0.95), names = FALSE)
      mc_lower[j] <- q[1]
      mc_upper[j] <- q[2]
    }
  }

  if (data_out) {
    return(data.frame(time = inp$years, type = inp$type, p_extant = p_mid,
                      lower = p_low, upper = p_upp,
                      mc_lower = mc_lower, mc_upper = mc_upper,
                      stringsAsFactors = FALSE))
  }

  n <- inp$n
  structure(
    list(p_extant = p_mid[n], lower = p_low[n], upper = p_upp[n],
         mc_lower = mc_lower[n], mc_upper = mc_upper[n],
         p_extinct = 1 - p_mid[n], end_year = inp$years[n],
         method = "Thompson et al. (2017) iterative model"),
    class = "ede_extant"
  )
}

#' @export
print.ede_extant <- function(x, ...) {
  cat(sprintf("<%s>\n", x$method))
  cat(sprintf("  probability extant at the end of %s: %s\n",
              format(x$end_year), format(signif(x$p_extant, 4))))
  cat(sprintf("  interval-arithmetic bounds: [%s, %s]\n",
              format(signif(x$lower, 4)), format(signif(x$upper, 4))))
  if (!is.na(x$mc_lower)) {
    cat(sprintf("  Monte Carlo 5%%-95%% range: [%s, %s]\n",
                format(signif(x$mc_lower, 4)), format(signif(x$mc_upper, 4))))
  }
  invisible(x)
}

#' Validate and expand Thompson et al. (2017, 2019) inputs to one row per year
#'
#' Returns the years, the type of each year and, for each input, vectors of
#' lower and upper bounds aligned with the years (record inputs are 0 in
#' survey years and vice versa).
#' @keywords internal
#' @noRd
thompson_inputs <- function(records, surveys, passive, p_initial,
                            start_year, end_year) {
  check_prob <- function(x, what) {
    if (!is.numeric(x) || anyNA(x) || any(x < 0 | x > 1)) {
      stop(sprintf("%s must contain values in [0, 1].", what), call. = FALSE)
    }
  }
  # bounds for a parameter from either `name` or `name_lower`/`name_upper`
  bounds <- function(x, name, what) {
    lo <- paste0(name, "_lower")
    hi <- paste0(name, "_upper")
    if (all(c(lo, hi) %in% names(x))) {
      out <- list(lower = unname(unlist(x[[lo]])), upper = unname(unlist(x[[hi]])))
    } else if (name %in% names(x)) {
      v <- unname(unlist(x[[name]]))
      out <- list(lower = v, upper = v)
    } else {
      stop(sprintf("%s must contain `%s` or `%s` and `%s`.", what, name, lo, hi),
           call. = FALSE)
    }
    check_prob(out$lower, sprintf("`%s` in %s", name, what))
    check_prob(out$upper, sprintf("`%s` in %s", name, what))
    if (any(out$lower > out$upper)) {
      stop(sprintf("lower bounds of `%s` in %s exceed upper bounds.", name, what),
           call. = FALSE)
    }
    out
  }
  check_years <- function(x, what) {
    if (is.null(x)) return(numeric(0))
    if (!is.data.frame(x) || !"year" %in% names(x)) {
      stop(sprintf("%s must be a data frame with a `year` column.", what),
           call. = FALSE)
    }
    y <- x$year
    if (!is.numeric(y) || anyNA(y) || any(y != round(y))) {
      stop(sprintf("`year` in %s must contain whole numbers.", what), call. = FALSE)
    }
    if (anyDuplicated(y)) {
      stop(sprintf("`year` in %s must not contain duplicates.", what), call. = FALSE)
    }
    y
  }

  ry <- check_years(records, "`records`")
  sy <- check_years(surveys, "`surveys`")
  if (length(ry) + length(sy) == 0L && (is.null(start_year) || is.null(end_year))) {
    stop("`start_year` and `end_year` are required when there are no records ",
         "or surveys.", call. = FALSE)
  }
  if (any(ry %in% sy)) {
    stop("a year cannot be both a record and an unsuccessful survey ",
         "(successful surveys are records).", call. = FALSE)
  }
  rb <- if (length(ry)) bounds(records, "pci", "`records`") else NULL
  sb <- if (length(sy)) {
    list(eps = bounds(surveys, "eps", "`surveys`"),
         pi = bounds(surveys, "pi", "`surveys`"),
         pr = bounds(surveys, "pr", "`surveys`"))
  } else NULL

  if (missing(passive) || is.null(passive)) {
    stop("`passive` must be supplied; use c(eps = 0, pi = 0, pr = 0) for ",
         "no passive surveillance.", call. = FALSE)
  }
  passive <- as.list(passive)
  pb <- list(eps = bounds(passive, "eps", "`passive`"),
             pi = bounds(passive, "pi", "`passive`"),
             pr = bounds(passive, "pr", "`passive`"))
  if (any(lengths(lapply(pb, `[[`, "lower")) != 1L)) {
    stop("`passive` must give one value (or one range) per input.", call. = FALSE)
  }

  if (missing(p_initial) || !length(p_initial) %in% 1:2) {
    stop("`p_initial` must be one probability or a c(lower, upper) range.",
         call. = FALSE)
  }
  check_prob(p_initial, "`p_initial`")
  if (p_initial[1] > p_initial[length(p_initial)]) {
    stop("`p_initial` must be given as c(lower, upper).", call. = FALSE)
  }
  p0 <- list(lower = p_initial[1], upper = p_initial[length(p_initial)])

  all_years <- c(ry, sy)
  check_year_arg <- function(x, what) {
    if (length(x) != 1L || !is.numeric(x) || is.na(x) || x != round(x)) {
      stop(sprintf("`%s` must be a whole number.", what), call. = FALSE)
    }
  }
  if (is.null(start_year)) start_year <- min(all_years)
  if (is.null(end_year)) end_year <- max(all_years)
  check_year_arg(start_year, "start_year")
  check_year_arg(end_year, "end_year")
  if (end_year < start_year) {
    stop("`end_year` must not be earlier than `start_year`.", call. = FALSE)
  }

  years <- seq(start_year, end_year)
  n <- length(years)
  type <- rep("passive", n)
  type[years %in% sy] <- "survey"
  type[years %in% ry] <- "record"

  fill <- function(idx_years, src, default) {
    lower <- upper <- rep(default, n)
    if (!is.null(src)) {
      m <- match(years, idx_years)
      ok <- !is.na(m)
      lower[ok] <- src$lower[m[ok]]
      upper[ok] <- src$upper[m[ok]]
    }
    list(lower = lower, upper = upper)
  }
  survey_input <- function(name) {
    x <- fill(sy, sb[[name]], NA_real_)
    gap <- type == "passive"
    x$lower[gap] <- pb[[name]]$lower
    x$upper[gap] <- pb[[name]]$upper
    x$lower[type == "record"] <- 0
    x$upper[type == "record"] <- 0
    x
  }
  pci <- fill(ry, rb, 0)

  list(years = years, type = type, n = n, p0 = p0, pci = pci,
       eps = survey_input("eps"), pi = survey_input("pi"),
       pr = survey_input("pr"))
}
