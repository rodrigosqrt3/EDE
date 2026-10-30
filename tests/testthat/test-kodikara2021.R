# Black-footed ferret: months since January 1972 (Kodikara 2021 data file
# kodikara2021_ferret_dates.xlsx, from Solow 1993; 28 sightings to
# December 1990; the authors' code uses T = 229 months).
k21_ferret <- sighting_data(data.frame(
  month = c(0, 6, 7, 8, 10, 17:22, 30, 31, 41, 44, 46, 53, 57, 58, 66, 90,
            117, 118, 122, 123, 127, 139, 151, 153),
  count = c(0, rep(1, 28))))
# Ivory-billed Woodpecker (authors' R_script_IBW.R): times are years since
# 1897; certain = physical evidence, uncertain = all other sightings.
k21_ibw_c <- sighting_data(data.frame(
  year = c(1897, 1898:1902, 1904:1910, 1913, 1914, 1917, 1924, 1925, 1932,
           1935, 1938, 1939),
  count = c(0, rep(1, 21))))
k21_ibw_u <- sighting_data(data.frame(
  year = c(1911, 1916, 1920, 1921, 1923, 1926, 1929:1931, 1933, 1934, 1936,
           1937, 1941:1944, 1946, 1948:1952, 1955, 1958, 1959, 1962,
           1966:1969, 1971:1974, 1976, 1981, 1982, 1985:1988, 1999,
           2004:2006),
  count = 1))
k21_cdf_at <- function(fit, t) {
  stats::approx(fit$density$time, fit$density$cdf, t, ties = "ordered")$y
}

test_that("black-footed ferret reproduces Table 1", {
  fit <- kodikara2021(k21_ferret, test_year = 229, n_grid = 100)
  yr <- 1972 + fit$extinction_time / 12
  # Table 1 / text: median March 1985, 95% HDI 1985 to March 1987 (JAGS)
  expect_equal(floor((yr[["median"]] - 1985) * 12), 2)   # March 1985
  expect_equal(round(yr[["lower"]]), 1985)
  expect_lt(abs(yr[["upper"]] - (1987 + 2.5 / 12)), 0.25)
  # Table 1: alpha_c median 0.71, 95% HDI (0.45, 1.01)
  expect_equal(round(unname(fit$shape["alpha_c", ]), 2), c(0.71, 0.45, 1.01),
               tolerance = 0.01)
  expect_lt(fit$posterior, 0.05)                        # extinct by 1990
  expect_s3_class(fit, "ede_bayes")
  expect_equal(fit$prior, exp(-0.005 * 229))
  expect_equal(fit$bayes_factor, (fit$posterior / (1 - fit$posterior)) /
                 (fit$prior / (1 - fit$prior)))
})

test_that("Ivory-billed Woodpecker reproduces Table 2, Figure 6 and Figure 7", {
  fit <- kodikara2021(k21_ibw_c, k21_ibw_u, test_year = 2010, n_grid = 100)
  # Table 2: 95% HDI 1939-1956, median 1950 (JAGS)
  expect_equal(round(unname(fit$extinction_time)), c(1950, 1939, 1956))
  # Figure 6: 95% HDIs alpha_c (0.37, 0.94), alpha_u1 (1.35, 3.46),
  # alpha_u2 (0.31, 1.92)
  hdi <- fit$shape[, c("lower", "upper")]
  expect_lt(max(abs(hdi - rbind(c(0.37, 0.94), c(1.35, 3.46), c(0.31, 1.92)))),
            0.04)
  # Figure 7 / text: P(tau_E <= 1948 | Dc, Du) = 0.38
  expect_lt(abs(k21_cdf_at(fit, 1948) - 0.38), 0.015)
  expect_lt(fit$posterior, 1e-6)
  # text: P(tau_E <= 1948 | Dc, alpha_c = 1) = 0.97
  hom <- kodikara2021(k21_ibw_c, test_year = 2010, homogeneous = TRUE,
                      n_grid = 50)
  expect_lt(abs(k21_cdf_at(hom, 1948) - 0.97), 0.01)
  expect_null(hom$shape)
  expect_match(hom$method, "homogeneous")
})

test_that("equations formulation differs from the published one", {
  fit <- kodikara2021(k21_ibw_c, k21_ibw_u, test_year = 2010,
                      formulation = "equations", n_grid = 60)
  expect_equal(fit$formulation, "equations")
  expect_lt(fit$extinction_time[["median"]], 1948)
  # eq. 9 leaves the likelihood constant beyond T: the HDI extends into the
  # exponential prior tail, computed analytically
  expect_gt(fit$extinction_time[["upper"]], 2010)
  expect_gt(fit$posterior, 0.2)
  expect_equal(max(fit$density$cdf), 1 - fit$posterior)
  # certain sightings only, equations
  fer <- kodikara2021(k21_ferret, test_year = 229, formulation = "equations",
                      n_grid = 60)
  expect_lt(fer$extinction_time[["median"]], 160)
})

test_that("upper incomplete gamma matches numerical integration", {
  ref <- function(s, x) {
    # int_x^Inf v^(s-1) e^-v dv with v = x e^t
    log(stats::integrate(function(t) exp(s * (log(x) + t) - x * exp(t)),
                         0, Inf, rel.tol = 1e-12)$value)
  }
  cases <- rbind(c(2.5, 0.3), c(-0.5, 0.2), c(-3.2, 0.05), c(-2, 0.5),
                 c(0, 0.3), c(-0.7, 3), c(-45, 0.5), c(0.4, 2))
  got <- kodikara2021_lgu(cases[, 1], log(cases[, 2]))
  want <- apply(cases, 1, function(r) ref(r[1], r[2]))
  expect_equal(got, want, tolerance = 1e-8)
  # underflowing and overflowing x
  expect_equal(kodikara2021_lgu(2, -700), lgamma(2))
  expect_equal(kodikara2021_lgu(-0.5, -700), -0.5 * -700 - log(0.5))
  expect_equal(kodikara2021_lgu(0, -700), log(700 - 0.5772156649015329))
  expect_equal(kodikara2021_lgu(-0.5, 800), -Inf)
})

test_that("block integrals match direct integration over (alpha, sigma)", {
  spec <- list(alpha = exp(seq(log(1e-4), log(5), length.out = 300)),
               smax = 30, homogeneous = FALSE)
  brute <- function(y, W, m) {
    K <- length(y)
    f_alpha <- function(a) vapply(a, function(al) {
      stats::integrate(function(s) {
        exp(K * log(al) - al * K * log(s) + (al - 1) * sum(log(y)) -
              m * W(al) * s^(-al))
      }, 0, spec$smax, rel.tol = 1e-10)$value / spec$smax
    }, numeric(1))
    log(stats::integrate(f_alpha, 0, 5, rel.tol = 1e-8,
                         subdivisions = 500L)$value / 5)
  }
  y <- c(1, 3, 4)
  got <- kodikara2021_block_total(y, function(a) a * log(6), 3, spec)$l
  expect_equal(got, brute(y, function(a) 6^a, 3), tolerance = 1e-4)
  # no sightings after tau (post-extinction block with K = 0)
  got0 <- kodikara2021_block_total(numeric(0),
                                   function(a) a * log(10) + log1p(-0.6^a),
                                   1, spec)$l
  # (flat integrand in alpha: trapezoid error on 300 log-spaced nodes ~1e-4)
  expect_equal(got0, brute(numeric(0), function(a) 10^a - 6^a, 1),
               tolerance = 1e-3)
  # zero exposure
  expect_equal(kodikara2021_block(1, 0, 0, -Inf, 30), 0)
  expect_equal(kodikara2021_block(1, 2, 0, -Inf, 30), -Inf)
  # homogeneous block: closed form with alpha = 1, K = 2, c = 6:
  # (1/30) int_0^30 s^-2 exp(-6/s) ds = exp(-6/30) / (6 * 30)
  hs <- modifyList(spec, list(homogeneous = TRUE))
  expect_equal(kodikara2021_block_total(c(1, 2), function(a) a * log(6), 1,
                                        hs)$l, -0.2 - log(180))
})

test_that("uncertain sightings only after the last certain one", {
  c_sd <- sighting_data(data.frame(t = c(0, 1, 2, 4, 5), n = c(0, 1, 2, 1, 1)))
  u_sd <- sighting_data(data.frame(t = c(7, 9, 12), n = 1))
  pub <- kodikara2021(c_sd, u_sd, test_year = 15, n_grid = 30,
                      sigma_max = 50, alpha_max = 20, rate = 0.05)
  expect_equal(nrow(pub$shape), 3)
  expect_true(all(diff(pub$density$cdf) >= -1e-12))
  expect_gt(pub$extinction_time[["lower"]], 5 - 1e-9)
  # zero-count uncertain record is the certain-only model
  none <- kodikara2021(c_sd, sighting_data(data.frame(t = 3, n = 0)),
                       test_year = 15, n_grid = 30)
  only <- kodikara2021(c_sd, test_year = 15, n_grid = 30)
  expect_equal(none$extinction_time, only$extinction_time)
  expect_equal(rownames(only$shape), "alpha_c")
  expect_output(print(pub), "Kodikara et al. \\(2021\\)")
})

test_that("input validation", {
  expect_error(kodikara2021(data.frame(t = 1, n = 1), test_year = 5),
               "sighting_data")
  expect_error(kodikara2021(k21_ferret, uncertain = 1, test_year = 300),
               "sighting_data")
  expect_error(kodikara2021(k21_ferret, test_year = 300, rate = 0), "rate")
  expect_error(kodikara2021(k21_ferret, test_year = 300, sigma_max = -1),
               "sigma_max")
  expect_error(kodikara2021(k21_ferret, test_year = 300, alpha_max = NA),
               "alpha_max")
  expect_error(kodikara2021(k21_ferret, test_year = 300, n_grid = "a"),
               "n_grid")
  expect_error(kodikara2021(k21_ferret, test_year = 300, level = 2), "level")
  expect_error(kodikara2021(k21_ferret, test_year = 300, homogeneous = NA),
               "homogeneous")
  expect_error(kodikara2021(k21_ibw_c, k21_ibw_u, test_year = 2010,
                            homogeneous = TRUE), "only available")
  expect_error(kodikara2021(sighting_data(data.frame(t = 1:3, n = 1)),
                            test_year = 5), "later than the origin")
  expect_error(kodikara2021(k21_ferret, test_year = 150), "later than")
  expect_error(kodikara2021(k21_ibw_c, k21_ibw_u, test_year = 2005),
               "not after")
  expect_error(kodikara2021(k21_ibw_c, sighting_data(data.frame(t = 1897, n = 1)),
                            test_year = 2010), "after the origin")
  expect_error(kodikara2021(k21_ferret, test_year = 300, formulation = "x"))
})
