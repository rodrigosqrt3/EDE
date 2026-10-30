#' Thompson et al. (2019) Bayesian updating of the probability of persistence
#'
#' Bayesian counterpart of [thompson2017()]. Each year's posterior is the
#' next year's prior (Thompson et al., 2019, equations 1-6), so the
#' probability that the taxon is extant after year `t` is
#' `P(X[t+1]) = 1 / (1 + B[t] * (1 / P(X[1]) - 1))` (equation 7), where
#' `P(X[1])` is the initial probability at the beginning of `start_year` and
#' `B[t]` is the cumulative Bayes factor for extinction, the product of the
#' yearly factors `b[j]` (equation 8). The yearly factors (equations 9-11) are
#'
#' * record year: `b = 1 - p(ci)`;
#' * unsuccessful dedicated survey year: `b = 1 / (1 - eps * p(r) * p(i))`;
#' * gap year (passive survey): the survey factor with `eps'`, `p'(r)` and
#'   `p'(i)`.
#'
#' Inputs use the same format as [thompson2017()]. Following section 2.1 and
#' the case study of the paper, ranges are replaced by their midpoints
#' (including a range for `p_initial`); the paper does not propagate the
#' ranges any further.
#'
#' The result is an `ede_bayes` object, as returned by [solow1993_bayes()].
#' Its `bayes_factor` is in favour of persistence, `1 / B[t]`, so that it has
#' the same orientation as in [solow1993_bayes()]; the paper's cumulative
#' Bayes factor for extinction is returned as `bayes_factor_extinction`.
#'
#' @inheritParams thompson2017
#' @param p_initial Probability `P(X[1])` that the taxon is extant at the
#'   beginning of `start_year`: one value or a `c(lower, upper)` range, whose
#'   midpoint is used.
#' @param data_out If `TRUE`, return the year-by-year Bayes factors and
#'   probabilities instead of the summary object.
#'
#' @return A list of class `ede_bayes` with components `bayes_factor`
#'   (cumulative Bayes factor in favour of persistence, `1 / B`),
#'   `bayes_factor_extinction` (`B`, equation 8), `posterior` (probability
#'   that the taxon is extant at the end of `end_year`, equation 7), `prior`
#'   (`P(X[1])`), `test_year` (`end_year`) and `method`; or (if
#'   `data_out = TRUE`) a data frame with one row per year and columns
#'   `time`, `type` (`"record"`, `"survey"` or `"passive"`), `b` (yearly
#'   Bayes factor for extinction, equation 11), `bayes_factor_extinction`
#'   (cumulative `B[t]`) and `posterior` (probability extant at the end of
#'   year `time`).
#'
#' @references
#' Thompson, C. J., Kodikara, S., Burgman, M. A., Demirhan, H., & Stone, L.
#'   (2019). Bayesian updating to estimate extinction from sequential
#'   observation data. Biological Conservation, 229, 26-29.
#'
#' @examples
#' # Alaotra Grebe from 1990 (Thompson et al., 2019, Table 1)
#' sv <- data.frame(year = c(1990, 1993, 1994, 1997, 1998, 1999),
#'                  eps_lower = 0.80, eps_upper = 0.95, pi_lower = 0.90,
#'                  pi_upper = 0.95, pr_lower = 0.70, pr_upper = 0.90)
#' pas <- c(eps_lower = 0, eps_upper = 0.05, pi_lower = 0.10, pi_upper = 0.65,
#'          pr_lower = 0.40, pr_upper = 0.60)
#' thompson2019(NULL, sv, pas, p_initial = 0.5, end_year = 1994)
#'
#' @export
thompson2019 <- function(records, surveys, passive, p_initial = 0.5,
                         start_year = NULL, end_year = NULL,
                         data_out = FALSE) {
  inp <- thompson_inputs(records, surveys, passive, p_initial,
                         start_year, end_year)
  mid <- function(x) (x$lower + x$upper) / 2
  p1 <- mid(inp$p0)
  if (p1 <= 0 || p1 >= 1) {
    stop("the (midpoint) initial probability must be in (0, 1).", call. = FALSE)
  }
  pci <- mid(inp$pci)
  eps_pr_pi <- mid(inp$eps) * mid(inp$pr) * mid(inp$pi)

  b <- ifelse(inp$type == "record", 1 - pci, 1 / (1 - eps_pr_pi))  # eq. 11
  B <- cumprod(b)                                                  # eq. 8
  if (anyNA(B)) {
    stop("a certain record (p(ci) = 1) and a certain unsuccessful survey ",
         "(eps * p(r) * p(i) = 1) cannot both occur.", call. = FALSE)
  }
  posterior <- 1 / (1 + B * (1 / p1 - 1))                          # eq. 7

  if (data_out) {
    return(data.frame(time = inp$years, type = inp$type, b = b,
                      bayes_factor_extinction = B, posterior = posterior,
                      stringsAsFactors = FALSE))
  }

  n <- inp$n
  structure(
    list(bayes_factor = 1 / B[n], bayes_factor_extinction = B[n],
         posterior = posterior[n], prior = p1, test_year = inp$years[n],
         method = "Thompson et al. (2019) Bayesian updating"),
    class = "ede_bayes"
  )
}
