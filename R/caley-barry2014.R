#' Caley & Barry (2014) Bayesian extinction probability from a detection record
#'
#' Posterior probability that a population is extant at the end of
#' `test_year`, from a record of detections (`y[t] = 1`) and non-detections
#' (`y[t] = 0`) in each period `t = 1, ..., T` (Caley & Barry, 2014). The
#' population is extant at the start of period 1. The latent variable `z` is
#' the period in which it goes extinct (`z > T` if it persists beyond the
#' record). If it goes extinct in period `z`, it cannot be detected from
#' period `z` on.
#'
#' **Constant model** (`model = "constant"`, section "Model specification -
#' constant population"). Each period, the population goes extinct with
#' probability `phi` and, if extant, is detected with probability `lambda`.
#' The priors are `phi ~ Beta(a_phi, b_phi)` and
#' `lambda ~ Beta(a_lambda, b_lambda)`. The paper's likelihood is
#' `[y | z] = prod_{t < z} lambda^y[t] (1 - lambda)^(1 - y[t])`, and
#' `[z] = (1 - phi)^(z - 1) phi`, or `(1 - phi)^T` when `z > T`. The paper
#' samples the posterior with a Gibbs sampler (Code S1). Here it is computed
#' exactly instead: given `z`, the integrals over `phi` and `lambda` are Beta
#' functions, so the posterior of `z` is a finite set of weights, and the
#' posteriors of `phi` and `lambda` are mixtures of Beta distributions.
#'
#' The published numbers come from Code S1. Its `lambda` update,
#' `Beta(a_lambda + S, b_lambda + z - S)`, counts the extinction period
#' itself as a non-detection, one period more than the likelihood above. The
#' `z` update has no such term. The sampler is still a valid Gibbs sampler:
#' the extra factor `(1 - lambda)` is the same for every `z`, so it targets
#' the model above with the prior `lambda ~ Beta(a_lambda, b_lambda + 1)`.
#' `formulation = "published"` (the default) computes that posterior and
#' reproduces the paper's results. `formulation = "equations"` uses the
#' likelihood as printed in the paper.
#'
#' **Non-constant model** (`model = "nonconstant"`, Code S2). The population
#' size is `N[t] = exp(r (t - 1))`, with `N[1] = 1` fixed because `N[1]` and
#' `delta` are not separately identifiable. Detection is
#' `lambda[t] = 1 - exp(-delta N[t])` (equation 3) and extinction is
#' `logit(phi[t]) = eps0 - eps1 log(N[t])` (equation 4). The priors are
#' uniform on `prior_growth`, `prior_delta`, `prior_eps0` and `prior_eps1`.
#' Each iteration of the Metropolis-within-Gibbs sampler of Code S2 is ported
#' as published:
#'
#' 1. impute `z` from its conditional distribution;
#' 2. take a log-normal random-walk step for `delta` (sd 0.5);
#' 3. take a joint step for `eps0` (normal, sd 1.5) and `eps1` (log-normal,
#'    sd 0.15);
#' 4. take a normal random-walk step for `r` (sd 0.05).
#'
#' The chain starts from `r = 0`, `delta = log(2)`, `eps0 = 0` and
#' `eps1 = 1`. As in Code S1, the `lambda` likelihood used in the `delta` and
#' `r` updates includes the extinction period as a non-detection, while the
#' `z` update does not. Because `lambda[t]` now varies with `t`, this factor
#' differs between values of `z`. **The full conditionals of Code S2 are
#' therefore not those of a single joint model**, and the sampler is
#' reproduced only for comparability with the paper. Code S2 uses
#' `r ~ U(-2.3, 0.69)`, whereas the paper's text gives `U(-1.6, 0.69)`; the
#' code's value is the default here. Code S2 also does not report an implied
#' prior probability of persistence that can be checked, so `prior` and
#' `bayes_factor` are `NA` for this model.
#'
#' **Input.** Periods are the integer time units from the earliest time in
#' `sd` to `test_year`. A period is a detection if its count is positive;
#' periods not listed in `sd` are non-detections. Unlike in [solow1993()],
#' the first period is an observation period even if it has a sighting. To
#' start observation before the first sighting, add a zero-count row (e.g.
#' the 1999 start of the paper's fox example).
#'
#' @param sd A [sighting_data] object with whole-number times.
#' @param test_year Last period of observation. Must be later than the last
#'   sighting. The posterior refers to persistence beyond the end of this
#'   period, given no detections after the last sighting.
#' @param model `"constant"` (exact computation) or `"nonconstant"` (MCMC).
#' @param formulation For the constant model: `"published"` (Code S1, which
#'   reproduces the paper) or `"equations"` (likelihood as printed).
#' @param prior_phi,prior_lambda Beta prior parameters `c(a, b)` for the
#'   constant model.
#' @param prior_growth,prior_delta,prior_eps0,prior_eps1 Uniform prior bounds
#'   for the non-constant model.
#' @param n_iter,burn_in Total MCMC iterations and initial iterations
#'   discarded (non-constant model). Use [set.seed()] for reproducibility.
#' @param data_out If `TRUE` (constant model only), return the posterior
#'   probability of persistence at the end of each period from the one after
#'   the last sighting to `test_year`, assuming no further detections.
#'
#' @return A list of class `ede_bayes` with components:
#'
#' * `posterior`: posterior probability that the population is extant at
#'   the end of `test_year`.
#' * `prior`: implied prior probability of that event,
#'   `B(a_phi, b_phi + T) / B(a_phi, b_phi)`, for the constant model; `NA`
#'   for the non-constant model.
#' * `bayes_factor`: posterior odds over prior odds of persistence; `NA` for
#'   the non-constant model.
#' * `test_year` and `method`.
#' * `extinction_time`: a data frame with the posterior probability that
#'   extinction occurred in each period after the last sighting.
#' * Posterior medians and 95% equal-tailed intervals: `phi` and `lambda`
#'   (constant model), or `growth`, `delta`, `eps0` and `eps1` (non-constant
#'   model).
#' * For the non-constant model, also the chain after burn-in (`samples`)
#'   and the acceptance rates (`acceptance`).
#'
#' If `data_out = TRUE`, a data frame with columns `time`, `bayes_factor`
#' and `posterior` is returned instead.
#'
#' @references
#' Caley, P., & Barry, S. C. (2014). Quantifying extinction probabilities
#'   from sighting records: inference and uncertainties. PLoS ONE, 9(4),
#'   e95857.
#'
#' @examples
#' # Tasmanian red fox carcasses (Caley & Barry, 2014, Table 1), 2001 start
#' fox <- sighting_data(data.frame(year = c(2001, 2003, 2005, 2006),
#'                                 count = 1))
#' caley_barry2014(fox, test_year = 2012)
#' # 1999 start
#' fox99 <- sighting_data(data.frame(year = c(1999, 2001, 2003, 2005, 2006),
#'                                   count = c(0, 1, 1, 1, 1)))
#' caley_barry2014(fox99, test_year = 2012)$posterior
#'
#' @export
caley_barry2014 <- function(sd, test_year, model = c("constant", "nonconstant"),
                            formulation = c("published", "equations"),
                            prior_phi = c(1.1, 1.1), prior_lambda = c(1.1, 1.1),
                            prior_growth = c(-2.3, 0.69),
                            prior_delta = c(0.01, 4.6),
                            prior_eps0 = c(-20, 20), prior_eps1 = c(0, 20),
                            n_iter = 1e5, burn_in = n_iter / 10,
                            data_out = FALSE) {
  stopifnot(inherits(sd, "sighting_data"))
  model <- match.arg(model)
  formulation <- match.arg(formulation)
  if (!any(sd$count > 0)) stop("no sightings with count > 0.", call. = FALSE)
  origin <- min(sd$time)
  if (any(sd$time != round(sd$time))) {
    stop("caley_barry2014() requires whole-number times (one row per period).",
         call. = FALSE)
  }
  last_sight <- max(sd$time[sd$count > 0])
  check_test_year(test_year, last_sight)
  if (test_year != round(test_year)) {
    stop("`test_year` must be a whole number.", call. = FALSE)
  }
  y <- as.integer(expand_record(sd, test_year)$count > 0)

  if (model == "constant") {
    check_pair <- function(x, what) {
      if (!is.numeric(x) || length(x) != 2L || anyNA(x) || any(x <= 0)) {
        stop(sprintf("`%s` must be two positive Beta parameters.", what),
             call. = FALSE)
      }
    }
    check_pair(prior_phi, "prior_phi")
    check_pair(prior_lambda, "prior_lambda")
    shift <- if (formulation == "published") 1 else 0
    if (data_out) {
      times <- candidate_times(last_sight, test_year)
      res <- lapply(times, function(t) {
        fit <- caley_barry_constant(y[seq_len(t - origin + 1)], prior_phi,
                                    prior_lambda, shift, summaries = FALSE)
        c(fit$bayes_factor, fit$posterior)
      })
      res <- do.call(rbind, res)
      return(data.frame(time = times, bayes_factor = res[, 1],
                        posterior = res[, 2]))
    }
    fit <- caley_barry_constant(y, prior_phi, prior_lambda, shift)
    fit$extinction_time$time <- origin + fit$extinction_time$time - 1
    method <- sprintf("Caley & Barry (2014) constant model (%s)", formulation)
    out <- c(list(bayes_factor = fit$bayes_factor, posterior = fit$posterior,
                  prior = fit$prior, test_year = test_year, method = method),
             fit[c("extinction_time", "phi", "lambda")])
    return(structure(out, class = "ede_bayes"))
  }

  if (data_out) {
    stop("`data_out` is only available for the constant model.", call. = FALSE)
  }
  check_bounds <- function(x, what, lower = -Inf) {
    if (!is.numeric(x) || length(x) != 2L || anyNA(x) || x[1] >= x[2] ||
        x[1] < lower) {
      stop(sprintf("`%s` must be c(lower, upper) with lower < upper.", what),
           call. = FALSE)
    }
  }
  check_bounds(prior_growth, "prior_growth")
  check_bounds(prior_delta, "prior_delta", lower = 0)
  check_bounds(prior_eps0, "prior_eps0")
  check_bounds(prior_eps1, "prior_eps1", lower = 0)
  if (length(n_iter) != 1L || !is.numeric(n_iter) || is.na(n_iter) ||
      length(burn_in) != 1L || !is.numeric(burn_in) || is.na(burn_in) ||
      burn_in < 0 || n_iter - burn_in < 1) {
    stop("`n_iter` must exceed `burn_in` (>= 0).", call. = FALSE)
  }
  init <- c(0, log(2), 0, 1)
  if (init[1] <= prior_growth[1] || init[1] >= prior_growth[2] ||
      init[2] <= prior_delta[1] || init[2] >= prior_delta[2] ||
      init[3] <= prior_eps0[1] || init[3] >= prior_eps0[2] ||
      init[4] <= prior_eps1[1] || init[4] >= prior_eps1[2]) {
    stop("the starting values (r = 0, delta = log(2), eps0 = 0, eps1 = 1) ",
         "must lie inside the prior bounds.", call. = FALSE)
  }
  chain <- caley_barry_nonconstant(y, round(n_iter), init, prior_growth,
                                   prior_delta, prior_eps0, prior_eps1)
  acceptance <- colMeans(chain$accept)
  samples <- chain$samples[-seq_len(round(burn_in)), , drop = FALSE]
  summ <- function(x) {
    stats::setNames(stats::quantile(x, c(0.5, 0.025, 0.975), names = FALSE),
                    c("median", "lower", "upper"))
  }
  te <- samples[, "z"]
  k <- last_sight - origin + 1
  zs <- seq(k + 1, length(y))
  extinction_time <- data.frame(
    time = origin + zs - 1,
    probability = vapply(zs, function(z) mean(!is.na(te) & te == z),
                         numeric(1))
  )
  samples <- as.data.frame(samples)
  samples$z <- origin + samples$z - 1
  names(samples)[names(samples) == "z"] <- "extinction_time"
  structure(
    list(bayes_factor = NA_real_, posterior = mean(is.na(te)),
         prior = NA_real_, test_year = test_year,
         method = "Caley & Barry (2014) non-constant model (Code S2 sampler)",
         extinction_time = extinction_time,
         growth = summ(samples$growth), delta = summ(samples$delta),
         eps0 = summ(samples$eps0), eps1 = summ(samples$eps1),
         acceptance = acceptance, samples = samples),
    class = "ede_bayes"
  )
}

#' Exact posterior of the Caley & Barry (2014) constant model
#'
#' `y` is the 0/1 record for periods 1..T; `shift = 1` adds one period to
#' the second Beta parameter of `lambda` (the Code S1 formulation).
#' @keywords internal
#' @noRd
caley_barry_constant <- function(y, prior_phi, prior_lambda, shift,
                                 summaries = TRUE) {
  T <- length(y)
  n <- sum(y)
  k <- max(which(y == 1))
  ap <- prior_phi[1]
  bp <- prior_phi[2]
  al <- prior_lambda[1]
  bl <- prior_lambda[2] + shift
  z <- if (k < T) seq(k + 1, T) else integer(0)
  # extinct in period z: phi (1 - phi)^(z - 1) lambda^n (1 - lambda)^(z - 1 - n)
  # persists beyond T:   (1 - phi)^T lambda^n (1 - lambda)^(T - n)
  log_w <- c(lbeta(ap + 1, bp + z - 1) + lbeta(al + n, bl + z - 1 - n),
             lbeta(ap, bp + T) + lbeta(al + n, bl + T - n))
  w <- exp(log_w - max(log_w))
  w <- w / sum(w)
  posterior <- w[length(w)]
  prior <- exp(lbeta(ap, bp + T) - lbeta(ap, bp))
  bayes_factor <- (posterior / (1 - posterior)) / (prior / (1 - prior))
  out <- list(posterior = posterior, prior = prior, bayes_factor = bayes_factor)
  if (!summaries) return(out)

  phi_a <- c(rep(ap + 1, length(z)), ap)
  phi_b <- c(bp + z - 1, bp + T)
  lam_a <- rep(al + n, length(z) + 1)
  lam_b <- c(bl + z - 1 - n, bl + T - n)
  mix_quantile <- function(p, a, b) {
    stats::uniroot(function(x) sum(w * stats::pbeta(x, a, b)) - p,
                   c(0, 1), tol = 1e-12)$root
  }
  summ <- function(a, b) {
    stats::setNames(vapply(c(0.5, 0.025, 0.975), mix_quantile, numeric(1),
                           a = a, b = b),
                    c("median", "lower", "upper"))
  }
  c(out, list(
    extinction_time = data.frame(time = z, probability = w[seq_along(z)]),
    phi = summ(phi_a, phi_b),
    lambda = summ(lam_a, lam_b)
  ))
}

#' Port of the Caley & Barry (2014) Code S2 Metropolis-within-Gibbs sampler
#' @keywords internal
#' @noRd
caley_barry_nonconstant <- function(y, n_iter, init, prior_growth, prior_delta,
                                    prior_eps0, prior_eps1) {
  obs <- length(y)
  final <- max(which(y == 1))
  tt <- seq_len(obs) - 1
  after <- seq(final + 1, obs)
  lam_f <- function(p) 1 - exp(-p[2] * exp(p[1] * tt))
  phi_f <- function(p) {
    lp <- p[3] - p[4] * p[1] * tt                 # eps0 - eps1 log(N[t])
    ifelse(lp > 50, 1, stats::plogis(lp))
  }
  # Code S2 lnL.lam: periods up to and including the imputed extinction period
  lnl_lam <- function(l, z) {
    m <- if (is.na(z)) obs else z
    yy <- y[seq_len(m)]
    ll <- l[seq_len(m)]
    sum(log(ll[yy != 0])) + sum(log(1 - ll[yy != 1]))
  }
  # Code S2 lnL.phi
  lnl_phi <- function(f, z) {
    if (is.na(z)) sum(log(1 - f)) else sum(log(1 - f[seq_len(z - 1)])) + log(f[z])
  }
  lunif <- function(x, b) stats::dunif(x, b[1], b[2], log = TRUE)

  p <- init
  samples <- matrix(NA_real_, n_iter, 5,
                    dimnames = list(NULL, c("growth", "delta", "eps0", "eps1", "z")))
  accept <- matrix(0, n_iter, 3,
                   dimnames = list(NULL, c("delta", "eps", "growth")))
  for (i in seq_len(n_iter)) {
    l <- lam_f(p)
    f <- phi_f(p)
    # impute the extinction period (Code S2 calc.pars / impute.TE)
    cs <- c(0, cumsum(log(1 - f[after]) + log(1 - l[after])))
    w <- c(exp(cs[seq_along(after)] + log(f[after])), exp(cs[length(after) + 1]))
    w[!is.finite(w) | w < 0] <- 0
    if (sum(w[seq_along(after)]) <= 0) w[length(w)] <- 1
    j <- sample.int(length(w), 1L, prob = w)
    z <- if (j <= length(after)) final + j else NA_integer_

    # delta: log-normal random walk
    q <- p
    q[2] <- stats::rlnorm(1, log(p[2]), 0.5)
    lr <- lnl_lam(lam_f(q), z) + lunif(q[2], prior_delta) -
      lnl_lam(l, z) - lunif(p[2], prior_delta) +
      stats::dlnorm(p[2], log(q[2]), 0.5, log = TRUE) -
      stats::dlnorm(q[2], log(p[2]), 0.5, log = TRUE)
    if (!is.na(lr) && stats::runif(1) < exp(lr)) {
      p <- q
      accept[i, 1] <- 1
    }
    # eps0 (normal) and eps1 (log-normal), jointly
    q <- p
    q[3:4] <- c(stats::rnorm(1, p[3], 1.5), stats::rlnorm(1, log(p[4]), 0.15))
    lr <- lnl_phi(phi_f(q), z) + lunif(q[3], prior_eps0) + lunif(q[4], prior_eps1) -
      lnl_phi(phi_f(p), z) - lunif(p[3], prior_eps0) - lunif(p[4], prior_eps1) +
      stats::dlnorm(p[4], log(q[4]), 0.15, log = TRUE) -
      stats::dlnorm(q[4], log(p[4]), 0.15, log = TRUE)
    if (!is.na(lr) && stats::runif(1) < exp(lr)) {
      p <- q
      accept[i, 2] <- 1
    }
    # growth rate: normal random walk
    q <- p
    q[1] <- stats::rnorm(1, p[1], 0.05)
    lr <- lnl_lam(lam_f(q), z) + lnl_phi(phi_f(q), z) + lunif(q[1], prior_growth) -
      lnl_lam(lam_f(p), z) - lnl_phi(phi_f(p), z) - lunif(p[1], prior_growth)
    if (!is.na(lr) && stats::runif(1) < exp(lr)) {
      p <- q
      accept[i, 3] <- 1
    }
    samples[i, ] <- c(p, z)
  }
  list(samples = samples, accept = accept)
}
