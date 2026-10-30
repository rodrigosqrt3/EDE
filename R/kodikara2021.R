#' Kodikara et al. (2021) extinction time with non-homogeneous Poisson sightings
#'
#' Posterior distribution of the extinction time `tau_E` when certain and
#' uncertain sightings follow independent non-homogeneous Poisson processes
#' with Weibull (power-law) rates, and `tau_E` is a change-point for the
#' uncertain sightings (Kodikara, Demirhan, Wang & Stone, 2021). Times are
#' measured from the observation origin, which is the earliest time in
#' `certain`. The observation period is `(0, T)`.
#'
#' * Certain sightings occur before `tau_E` with rate
#'   `lambda_c(t) = (alpha_c / sigma_c) (t / sigma_c)^(alpha_c - 1)` and mean
#'   function `m_c(t) = (t / sigma_c)^alpha_c` (equations 3 and 5).
#' * Uncertain sightings occur with rate `lambda_u1` (parameters `alpha_u1`,
#'   `sigma_u1`) before `tau_E` and `lambda_u2` (`alpha_u2`, `sigma_u2`)
#'   after it (equations 4 and 6).
#' * The likelihood is equation 9: with `Kc` certain sightings at times
#'   `y_ci`, and `Ku(tau_E)` of the `Ku` uncertain sightings before `tau_E`,
#'   `prod_i lambda_c(y_ci) exp(-m_c(tau_E))`
#'   `x prod_(j <= Ku(tau_E)) lambda_u1(y_uj) exp(-m_u1(tau_E))`
#'   `x prod_(j > Ku(tau_E)) lambda_u2(y_uj) exp(-(m_u2(T) - m_u2(tau_E)))`;
#'   without uncertain sightings only the first factor is used
#'   (equation 10). The likelihood is zero for `tau_E` before the last
#'   certain sighting, and is evaluated at `tau_E = T` for `tau_E >= T`.
#' * Priors: `Uniform(0, alpha_max)` for the shape parameters,
#'   `Uniform(0, sigma_max)` for the scale parameters (both 1000 in the
#'   paper) and an exponential prior with rate `rate` (0.005 in the paper)
#'   for `tau_E` (equation 11). `homogeneous = TRUE` fixes all shape
#'   parameters at 1 (constant rates; section 3.2 and Figure 7).
#'
#' **Formulation.** The paper samples the posterior with JAGS. In the
#' authors' code (`R_script_Ferret.R`, `R_script_IBW.R`), the likelihood
#' is a product over sightings in which each sighting's factor contains
#' the exponential term of its block. The term `exp(-m_c(tau_E))` is
#' therefore raised to the power `Kc`, `exp(-m_u1(tau_E))` to the power
#' `Ku(tau_E)` and `exp(-(m_u2(T) - m_u2(tau_E)))` to the power
#' `Ku - Ku(tau_E)`. The code also uses `m_u1(tau_E)` without truncating
#' `tau_E` at `T`. `formulation = "published"` (the default) uses this
#' likelihood and reproduces the published results.
#' `formulation = "equations"` uses equation 9 as printed. Its results
#' differ from the published ones (for the Ivory-billed Woodpecker, the
#' posterior median of `tau_E` is 1946 instead of 1950), so it is not
#' verified against the paper.
#'
#' **Computation (package construct).** This function does not use MCMC.
#' Given `tau_E`, the likelihood factorises into the certain, pre- and
#' post-extinction uncertain blocks. In each block, the integral over the
#' scale `sigma` is an upper incomplete gamma function,
#' `int_0^sigma_max sigma^(-alpha K) exp(-c sigma^(-alpha)) d sigma =
#' c^(1/alpha - K) Gamma(K - 1/alpha, c sigma_max^(-alpha)) / alpha`.
#' The shape `alpha` is integrated numerically on a log-spaced grid. The
#' posterior density of `tau_E` is evaluated on a grid between the
#' discontinuities at the uncertain sighting times, and integrated with
#' the trapezoidal rule. The density beyond `T` is integrated exactly for
#' the equations (constant likelihood), and on a grid for the published
#' formulation. The median and the HDI (the shortest interval with
#' posterior probability `level`) are as reported in Tables 1 and 2 of the
#' paper. The marginal posteriors of the shape parameters are summarised in
#' the same way.
#'
#' **Input.** Times are in the units of `certain` and `uncertain`. The
#' priors depend on these units: the paper uses months for the
#' black-footed ferret and years for the Ivory-billed Woodpecker. A count
#' greater than one is that number of sightings at the same time. Every
#' sighting must be later than the origin. Add a zero-count row in
#' `certain` at the start of observation.
#'
#' @param certain A [sighting_data] object with the certain sightings. Its
#'   earliest time is the origin; at least one certain sighting must be
#'   later than it.
#' @param uncertain A [sighting_data] object with the uncertain sightings,
#'   or `NULL` (certain sightings only, equation 10). Zero-count rows are
#'   ignored.
#' @param test_year End `T` of the observation period. Must be later than
#'   the last certain sighting and not earlier than any uncertain sighting.
#' @param formulation `"published"` (the authors' code; reproduces the paper)
#'   or `"equations"` (equation 9 as printed).
#' @param homogeneous If `TRUE`, all shape parameters are fixed at 1. Only
#'   available without uncertain sightings: the paper's value for the
#'   homogeneous model with uncertain sightings (Figure 7) is not
#'   reproduced.
#' @param rate Rate of the exponential prior of `tau_E`, per time unit.
#' @param sigma_max,alpha_max Upper bounds of the uniform priors of the
#'   scale and shape parameters.
#' @param level Probability content of the HDIs.
#' @param n_grid Approximate number of grid points for `tau_E` in
#'   `(last certain sighting, T)`.
#'
#' @return A list of class `ede_bayes` with components:
#'
#' * `posterior`: posterior probability `P(tau_E > T)` of persistence beyond
#'   the observation period.
#' * `prior`: the prior probability of that event, `exp(-rate (T - origin))`.
#' * `bayes_factor`: posterior odds over prior odds of persistence.
#' * `test_year`, `method`, `formulation` and `level`.
#' * `extinction_time`: posterior median and HDI of `tau_E`, as times
#'   (`c(median, lower, upper)`).
#' * `shape`: posterior medians and HDIs of `alpha_c`, `alpha_u1` and
#'   `alpha_u2` (a matrix with one row for each estimated shape).
#' * `density`: data frame with the posterior density (`density`) and
#'   distribution function (`cdf`) of `tau_E` on the grid (`time`). Use it
#'   for `P(tau_E <= t)`.
#'
#' @references
#' Kodikara, S., Demirhan, H., Wang, Y. & Stone, L. (2021). Inferring
#'   extinction date of a species using non-homogeneous Poisson processes
#'   with a change-point. Methods in Ecology and Evolution, 12(3), 530-538.
#'   doi:10.1111/2041-210X.13542
#'
#' @examples
#' # Black-footed ferret, months since January 1972 (Kodikara et al., 2021,
#' # Table 1; dates from Solow, 1993)
#' months <- (ferret$year - 1972) * 12 + ferret$month
#' sd <- sighting_data(data.frame(month = c(0, months),
#'                                count = c(0, ferret$count)))
#' fit <- kodikara2021(sd, test_year = 229, n_grid = 200)
#' 1972 + fit$extinction_time / 12
#' fit$shape
#'
#' @export
kodikara2021 <- function(certain, uncertain = NULL, test_year,
                         formulation = c("published", "equations"),
                         homogeneous = FALSE, rate = 0.005, sigma_max = 1000,
                         alpha_max = 1000, level = 0.95, n_grid = 200) {
  formulation <- match.arg(formulation)
  if (!inherits(certain, "sighting_data")) {
    stop("`certain` must be a sighting_data object.", call. = FALSE)
  }
  if (!is.null(uncertain) && !inherits(uncertain, "sighting_data")) {
    stop("`uncertain` must be a sighting_data object or NULL.", call. = FALSE)
  }
  positive <- function(x, what) {
    if (length(x) != 1L || !is.numeric(x) || is.na(x) || !is.finite(x) ||
        x <= 0) {
      stop(sprintf("`%s` must be a positive number.", what), call. = FALSE)
    }
  }
  positive(rate, "rate")
  positive(sigma_max, "sigma_max")
  positive(alpha_max, "alpha_max")
  positive(n_grid, "n_grid")
  if (length(level) != 1L || !is.numeric(level) || is.na(level) ||
      level <= 0 || level >= 1) {
    stop("`level` must be in (0, 1).", call. = FALSE)
  }
  if (!is.logical(homogeneous) || length(homogeneous) != 1L ||
      is.na(homogeneous)) {
    stop("`homogeneous` must be TRUE or FALSE.", call. = FALSE)
  }
  if (homogeneous && !is.null(uncertain)) {
    stop("`homogeneous = TRUE` is only available without uncertain ",
         "sightings (not verified against the paper otherwise).",
         call. = FALSE)
  }
  origin <- min(certain$time)
  keep <- certain$count > 0
  yc <- rep(certain$time[keep], certain$count[keep]) - origin
  if (!length(yc) || any(yc <= 0)) {
    stop("certain sightings must be later than the origin (the earliest ",
         "time in `certain`); add a zero-count row at the start of ",
         "observation.", call. = FALSE)
  }
  last_certain <- max(yc) + origin
  check_test_year(test_year, last_certain)
  yu <- NULL
  if (!is.null(uncertain)) {
    keep <- uncertain$count > 0
    yu <- sort(rep(uncertain$time[keep], uncertain$count[keep]) - origin)
    if (!length(yu)) {
      yu <- NULL
    } else if (any(yu <= 0) || any(yu > test_year - origin)) {
      stop("uncertain sightings must lie after the origin (the earliest ",
           "time in `certain`) and not after `test_year`.", call. = FALSE)
    }
  }

  spec <- list(yc = yc, yu = yu, T = test_year - origin, tL = max(yc),
               published = formulation == "published",
               homogeneous = homogeneous, rate = rate, smax = sigma_max,
               alpha = exp(seq(log(1e-4), log(alpha_max), length.out = 300)))
  post <- kodikara2021_posterior(spec, n_grid)

  # summaries of tau_E
  grid <- post$grid
  cdf <- post$cdf
  keep <- !duplicated(cdf)
  top <- max(cdf)
  inv <- function(p) {
    if (p <= top) return(stats::approx(cdf[keep], grid[keep], p)$y)
    # beyond the grid: exponential prior tail of the equations formulation
    if (spec$published || p >= 1) return(Inf)
    spec$T + log((1 - top) / (1 - p)) / rate
  }
  best <- c(NA_real_, NA_real_)
  width <- Inf
  for (i in seq_along(grid)) {
    up <- inv(cdf[i] + level)
    if (!is.finite(up)) break
    if (up - grid[i] < width) {
      width <- up - grid[i]
      best <- c(grid[i], up)
    }
  }
  ext <- origin + c(median = inv(0.5), lower = best[1], upper = best[2])

  pers <- post$persist
  prior <- exp(-rate * spec$T)
  bf <- (pers / (1 - pers)) / (prior / (1 - prior))
  shape <- NULL
  if (!homogeneous) {
    shape <- t(vapply(post$alpha_density, kodikara2021_shape_summary,
                      numeric(3), alpha = spec$alpha, level = level))
    rownames(shape) <- c("alpha_c", "alpha_u1", "alpha_u2")[seq_len(nrow(shape))]
  }
  structure(
    list(bayes_factor = bf, posterior = pers, prior = prior,
         test_year = test_year,
         method = sprintf("Kodikara et al. (2021) %s Poisson model (%s)",
                          if (homogeneous) "homogeneous" else "Weibull-rate",
                          formulation),
         formulation = formulation, level = level,
         extinction_time = ext, shape = shape,
         density = data.frame(time = origin + grid, density = post$density,
                              cdf = cdf)),
    class = "ede_bayes"
  )
}

#' log of the upper incomplete gamma function Gamma(s, x), any real s, given
#' log(x) (vectorised). s > 0 uses pgamma; s <= 0 a continued fraction
#' (x - s >= 1) or one downward recurrence step from s + 1 (s in (-1, 0]).
#' @keywords internal
#' @noRd
kodikara2021_lgu <- function(s, lx) {
  n <- max(length(s), length(lx))
  s <- rep_len(s, n)
  lx <- rep_len(lx, n)
  x <- exp(lx)
  out <- numeric(n)
  tiny <- lx < -600                               # x underflows
  pos <- s > 0
  i <- pos & !tiny
  out[i] <- lgamma(s[i]) +
    stats::pgamma(x[i], s[i], lower.tail = FALSE, log.p = TRUE)
  i <- pos & tiny
  out[i] <- lgamma(s[i])
  # x -> 0, s < 0: Gamma(s, x) ~ x^s / (-s); s = 0: E1(x) ~ -log(x) - gamma
  i <- !pos & tiny & s < 0
  out[i] <- s[i] * lx[i] - log(-s[i])
  i <- !pos & tiny & s == 0
  out[i] <- log(-lx[i] - 0.5772156649015329)
  big <- !pos & !tiny & !is.finite(x)
  out[big] <- -Inf
  neg <- which(!pos & !tiny & is.finite(x))
  if (!length(neg)) return(out)
  sn <- s[neg]
  xn <- x[neg]
  cf <- xn >= 1 | sn <= -30
  if (any(cf)) {
    # modified Lentz continued fraction (Legendre)
    ss <- sn[cf]
    xx <- xn[cf]
    fpmin <- 1e-300
    b <- xx + 1 - ss
    cc <- rep(1 / fpmin, length(b))
    d <- 1 / b
    h <- d
    for (k in seq_len(5000)) {
      an <- -k * (k - ss)
      b <- b + 2
      d <- an * d + b
      d[abs(d) < fpmin] <- fpmin
      cc <- b + an / cc
      cc[abs(cc) < fpmin] <- fpmin
      d <- 1 / d
      del <- d * cc
      h <- h * del
      if (all(abs(del - 1) < 1e-15)) break
    }
    out[neg[cf]] <- -xx + ss * log(xx) + log(h)
  }
  rest <- which(!cf)
  if (length(rest)) {
    # -30 < s <= 0 and x < 1: start at a0 = s + n in (0, 1] (or at
    # Gamma(0, x) = E1(x) when a0 is 1) and recur downwards,
    # Gamma(a - 1, x) = (x^(a - 1) exp(-x) - Gamma(a, x)) / (1 - a),
    # which is stable for x < 1.
    ss <- sn[rest]
    xx <- xn[rest]
    n <- floor(-ss) + 1
    a <- ss + n
    e1 <- abs(a - 1) < 1e-8
    k <- seq_len(60)
    E1 <- function(x) {
      -0.5772156649015329 - log(x) -
        vapply(x, function(z) sum((-z)^k / (k * factorial(k))), numeric(1))
    }
    g <- numeric(length(ss))
    g[!e1] <- exp(lgamma(a[!e1]) + stats::pgamma(xx[!e1], a[!e1],
                                                lower.tail = FALSE,
                                                log.p = TRUE))
    g[e1] <- E1(xx[e1])
    a[e1] <- 0
    n[e1] <- n[e1] - 1
    for (step in seq_len(max(n))) {
      go <- n >= step
      g[go] <- (xx[go]^(a[go] - 1) * exp(-xx[go]) - g[go]) / (1 - a[go])
      a[go] <- a[go] - 1
    }
    out[neg[rest]] <- log(g)
  }
  out
}

#' log of (1 / sigma_max) int_0^sigma_max alpha^K sigma^(-alpha K)
#'   prod(y^(alpha - 1)) exp(-c sigma^(-alpha)) d sigma, with log(c) = lc
#'   (vectors over alpha); S = sum(log(y)).
#' @keywords internal
#' @noRd
kodikara2021_block <- function(alpha, K, S, lc, smax) {
  out <- numeric(length(alpha))
  zero <- lc == -Inf                              # no exposure: only K = 0
  out[zero] <- if (K == 0) 0 else -Inf
  a <- alpha[!zero]
  l <- lc[!zero]
  out[!zero] <- (K - 1) * log(a) + (a - 1) * S + (1 / a - K) * l +
    kodikara2021_lgu(K - 1 / a, l - a * log(smax)) - log(smax)
  out
}

#' Log likelihood of one block integrated over sigma and alpha. logW(alpha)
#' is the log exposure (log of t^alpha or T^alpha - t^alpha), m the power of
#' the exponential term. Returns the total `l` and the log integrand `g` on
#' the alpha grid (for the marginal posterior of alpha).
#' @keywords internal
#' @noRd
kodikara2021_block_total <- function(y, logW, m, spec) {
  K <- length(y)
  S <- sum(log(y))
  alpha <- spec$alpha
  if (K == 0 && m == 0) return(list(l = 0, g = rep(0, length(alpha))))
  if (spec$homogeneous) {
    return(list(l = kodikara2021_block(1, K, S, log(m) + logW(1), spec$smax),
                g = NULL))
  }
  g <- kodikara2021_block(alpha, K, S, log(m) + logW(alpha), spec$smax)
  u <- log(alpha)
  wq <- alpha * (c(diff(u), 0) + c(0, diff(u))) / 2  # trapezoid in log(alpha)
  mx <- max(g)
  list(l = mx + log(sum(wq * exp(g - mx))) - log(max(alpha)), g = g)
}

#' Log posterior density (unnormalised) of tau_E and the block integrands.
#' `n_before` is the number of uncertain sightings before tau_E (fixed
#' within a segment between uncertain sighting times).
#' @keywords internal
#' @noRd
kodikara2021_logpost <- function(tau, n_before, spec) {
  T <- spec$T
  tc <- min(tau, T)
  Kc <- length(spec$yc)
  B <- list(c = kodikara2021_block_total(
    spec$yc, function(a) a * log(tc), if (spec$published) Kc else 1, spec))
  if (!is.null(spec$yu)) {
    pre <- spec$yu[seq_len(n_before)]
    post <- spec$yu[setdiff(seq_along(spec$yu), seq_len(n_before))]
    t1 <- if (spec$published) tau else tc        # the code does not truncate
    after <- function(a) a * log(T) + log1p(-(tc / T)^a)
    if (spec$published) {
      B$u1 <- kodikara2021_block_total(pre, function(a) a * log(t1),
                                       length(pre), spec)
      B$u2 <- kodikara2021_block_total(post, after, length(post), spec)
    } else {
      B$u1 <- kodikara2021_block_total(pre, function(a) a * log(t1), 1, spec)
      B$u2 <- if (tc >= T) {
        list(l = 0, g = rep(0, length(spec$alpha)))
      } else {
        kodikara2021_block_total(post, after, 1, spec)
      }
    }
  }
  lp <- sum(vapply(B, function(b) b$l, numeric(1))) + log(spec$rate) -
    spec$rate * tau
  list(lp = lp, B = B)
}

#' Posterior of tau_E on a grid: density, cdf, persistence probability and
#' the marginal posterior densities of the shape parameters.
#' @keywords internal
#' @noRd
kodikara2021_posterior <- function(spec, n_grid) {
  T <- spec$T
  tL <- spec$tL
  yu <- spec$yu
  breaks <- sort(unique(c(tL, yu[yu > tL & yu < T], T)))
  step <- (T - tL) / n_grid
  nodes <- list()
  for (i in seq_len(length(breaks) - 1L)) {
    a <- breaks[i]
    b <- breaks[i + 1L]
    n <- max(2L, ceiling((b - a) / step) + 1L)
    nb <- sum(yu <= (a + b) / 2)                  # constant on the segment
    nodes[[i]] <- cbind(seq(a, b, length.out = n), nb, i)
  }
  if (spec$published) {
    # beyond T the published likelihood still varies; integrate on a
    # log-spaced grid up to 50 prior e-folding times
    ext <- T + exp(seq(log(step / 10), log(50 / spec$rate), length.out = 200))
    nodes[[length(nodes) + 1L]] <- cbind(c(T, ext), length(yu), length(nodes) + 1L)
  }
  nodes <- do.call(rbind, nodes)
  res <- lapply(seq_len(nrow(nodes)), function(i) {
    kodikara2021_logpost(nodes[i, 1], nodes[i, 2], spec)
  })
  lp <- vapply(res, function(r) r$lp, numeric(1))
  mx <- max(lp)
  d <- exp(lp - mx)
  grid <- nodes[, 1]
  seg <- nodes[, 3]
  same <- c(seg[-1] == seg[-length(seg)], FALSE)
  h <- c(diff(grid), 0) * same
  piece <- h * (d + c(d[-1], 0)) / 2              # trapezoid within segments
  tw <- (h + c(0, h[-length(h)])) / 2             # node weights
  cm <- c(0, cumsum(piece[-length(piece)]))
  total <- sum(piece)
  if (!spec$published) {
    # constant likelihood for tau_E >= T: exact exponential tail
    total <- total + d[length(d)] / spec$rate
  }
  in_obs <- grid <= T
  alpha_density <- NULL
  if (!spec$homogeneous) {
    nm <- names(res[[1]]$B)
    alpha_density <- lapply(nm, function(b) {
      dens <- numeric(length(spec$alpha))
      for (i in seq_along(res)) {
        bl <- res[[i]]$B[[b]]
        if (tw[i] > 0) dens <- dens + tw[i] * exp(lp[i] - mx - bl$l + bl$g)
      }
      if (!spec$published) {
        bl <- res[[length(res)]]$B[[b]]
        dens <- dens + exp(lp[length(lp)] - mx - bl$l + bl$g) / spec$rate
      }
      dens
    })
    names(alpha_density) <- nm
  }
  list(grid = grid, density = d / total, cdf = cm / total,
       persist = 1 - max(cm[in_obs]) / total, alpha_density = alpha_density)
}

#' Median and HDI of a shape parameter from its unnormalised density on the
#' (log-spaced) alpha grid.
#' @keywords internal
#' @noRd
kodikara2021_shape_summary <- function(dens, alpha, level) {
  hi <- alpha[max(which(dens >= max(dens) * 1e-12))]
  x <- seq(0, hi, length.out = 20001)[-1]
  dx <- stats::approx(alpha, dens, x, rule = 2)$y
  dx <- dx / sum(dx)
  o <- order(dx, decreasing = TRUE)
  inside <- o[seq_len(which(cumsum(dx[o]) >= level)[1])]
  c(median = x[which(cumsum(dx) >= 0.5)[1]], lower = x[min(inside)],
    upper = x[max(inside)])
}
