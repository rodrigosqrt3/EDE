ferret_record <- function() {
  # Solow (1993b), Table 2: months with January 1972 = 1, observation from 0
  year <- c(rep(1972, 4), rep(1973, 6), rep(1974, 2), rep(1975, 3),
            rep(1976, 3), 1977, 1979, 1981, 1981, rep(1982, 3), 1983, 1984, 1984)
  month <- c(6, 7, 8, 10, 5:10, 6, 7, 5, 8, 10, 5, 9, 10, 6, 6, 9, 10, 2, 3, 7, 7, 7, 9)
  time <- (year - 1972) * 12 + month
  sighting_data(data.frame(time = c(0, time), count = c(0, rep(1, length(time)))))
}

test_that("solow1993b reproduces published value for Black-footed ferret (Solow 1993b)", {
  p <- EDE:::solow1993b_fisher_series(153, 1714, 28) / EDE:::solow1993b_fisher_series(229, 1714, 28)
  expect_equal(round(p, 3), 0.050)

  d <- ferret_record()
  expect_equal(sum(d$count), 28)
  expect_equal(sum(d$time * d$count), 1714)
  curve <- solow1993b(d, test_year = 229, data_out = TRUE)
  expect_equal(round(curve$chance[curve$time == 229], 3), 0.050)
})

test_that("solow1993b returns the inverted-test upper bound", {
  d <- ferret_record()
  res <- solow1993b(d, test_year = 600)
  expect_true(res$estimate > 153 && res$estimate < res$upper)
  expect_equal(res$interval_type, "one-sided")
  p_at_bound <- EDE:::solow1993b_fisher_series(153, 1714, 28) /
    EDE:::solow1993b_fisher_series(res$upper, 1714, 28)
  expect_lt(abs(p_at_bound - 0.05), 1e-6)
  # the paper finds p = 0.050 at T = 229
  expect_lt(abs(res$upper - 229), 1)
})

test_that("solow1993b reports a missing bound when F(t[n]) exceeds alpha", {
  d <- ferret_record()
  floor_p <- EDE:::solow1993b_fisher_series(153, 1714, 28)
  expect_warning(res <- solow1993b(d, alpha = floor_p / 2, test_year = 1e5),
                 "never falls")
  expect_true(is.na(res$upper))
})

test_that("solow1993b validates alpha, test_year and minimum n", {
  d <- sighting_data(data.frame(year = c(1900, 1910), sightings = c(1, 1)))
  expect_error(solow1993b(d, alpha = 0, test_year = 2000), "must be in \\(0, 1\\)")
  expect_error(solow1993b(d, alpha = 1, test_year = 2000), "must be in \\(0, 1\\)")
  expect_error(solow1993b(d, 0.05), "test_year")
  expect_error(solow1993b(d, 0.05, test_year = "2000"), "test_year")
  expect_error(solow1993b(d, 0.05, test_year = 1910), "later than")

  d1 <- sighting_data(data.frame(year = 1900, sightings = 1))
  expect_error(solow1993b(d1, 0.05, test_year = 2000), "at least 2")
})

test_that("solow1993b data_out = TRUE returns full curve", {
  d <- sighting_data(data.frame(
    years = c(1907, 1910, 1915, 1916, 1920, 1925, 1930, 1931),
    sightings = c(1, 1, 3, 4, 3, 1, 2, 1)
  ))
  curve <- solow1993b(d, alpha = 0.05, test_year = 1950, data_out = TRUE)
  expect_named(curve, c("time", "chance"))
  expect_equal(nrow(curve), 19L)
})

test_that("solow1993b_fisher_series handles edge cases", {
  expect_equal(EDE:::solow1993b_fisher_series(0, 100, 5), 0)
  expect_equal(EDE:::solow1993b_fisher_series(10, 0, 5), 0)
  expect_equal(EDE:::solow1993b_fisher_series(200, 100, 5), 1)

  d_single <- data.frame(time = 1900, count = 1)
  expect_equal(EDE:::solow1993b_chance(d_single), 1.0)

  d_one_used <- data.frame(time = 0:2, count = c(0, 1, 0))
  expect_equal(EDE:::solow1993b_chance(d_one_used), 1.0)

  d_zero_denominator <- data.frame(time = 0:1, count = c(0, 2))
  expect_equal(EDE:::solow1993b_chance(d_zero_denominator), 1.0)
})

test_that("solow1993b warns when alpha is not reached before test_year", {
  d <- sighting_data(data.frame(
    years = c(1907, 1910, 1915, 1916, 1920, 1925, 1930, 1931),
    sightings = c(1, 1, 3, 4, 3, 1, 2, 1)
  ))
  expect_warning(
    res_na <- solow1993b(d, alpha = 1e-10, test_year = 1932),
    "never falls"
  )
  expect_true(is.na(res_na$upper))
})

test_that("solow1993b reproduces the Snakeblenny example (Solow & Helser 2000)", {
  # Solow & Helser (2000), p. 5: observation 1969-1992, 11 sightings;
  # declining-population test p ~ 0.348 (constant-rate test p ~ 0.036)
  d <- sighting_data(data.frame(
    time = c(1969, 1970:1975, 1977, 1978, 1980, 1981, 1986),
    count = c(0, rep(1, 11))
  ))
  curve <- solow1993b(d, test_year = 1992, data_out = TRUE)
  expect_equal(round(curve$chance[curve$time == 1992], 3), 0.348)
  curve_a <- solow1993(d, test_year = 1992, data_out = TRUE)
  expect_equal(round(curve_a$chance[curve_a$time == 1992], 3), 0.036)
})

test_that("solow1993b estimate is Solow (2005) equation 8", {
  # direct evaluation of the alternating sums of eq. 8, fine for small cases
  eq8 <- function(tn, s, n) {
    i <- 0:floor(s / tn)
    j <- 0:(floor(s / tn) - 1)
    tn + sum((-1)^i * choose(n, i) * (s - i * tn)^(n - 1)) /
      (n * (n - 1) * sum((-1)^j * choose(n - 1, j) * (s - (j + 1) * tn)^(n - 2)))
  }
  for (x in list(c(153, 1714, 28), c(10, 25, 3), c(7, 40, 9), c(100, 1500, 60))) {
    expect_lt(abs(EDE:::solow1993b_mvue(x[1], x[2], x[3]) - eq8(x[1], x[2], x[3])), 1e-6)
  }
  # n = 2: eq. 8 reduces to t[n] + (2 t[n] - s) / 2
  expect_equal(EDE:::solow1993b_mvue(10, 15, 2), 12.5)

  # ferret record: estimate in months from January 1972 (exact value from
  # rational arithmetic: 164.302188809)
  res <- solow1993b(ferret_record(), test_year = 600)
  expect_lt(abs(res$estimate - 164.302188809), 1e-6)
})

test_that("solow1993b estimate stays accurate where the series cancels", {
  # exact values from rational arithmetic (gmp) of eq. 8; the alternating
  # series in double precision gives 592.9 and 50.29 for these cases
  expect_lt(abs(EDE:::solow1993b_mvue(100, 5900, 60) - 100.028248588), 1e-6)
  expect_lt(abs(EDE:::solow1993b_mvue(50, 4000, 200) - 50.4954943073), 1e-6)
  expect_lt(abs(EDE:::solow1993b_mvue(100, 3000, 60) - 101.694915254), 1e-6)
})

test_that("solow1993b estimate is unbiased under the truncated exponential model", {
  set.seed(42)
  b <- 1 / 60; te <- 150; n <- 15
  est <- replicate(4000, {
    x <- -log(1 - stats::runif(n) * (1 - exp(-b * te))) / b
    EDE:::solow1993b_mvue(max(x), sum(x), n)
  })
  expect_lt(abs(mean(est) - te), 4 * stats::sd(est) / sqrt(length(est)))
})

test_that("solow1993b_log_fisher matches Fisher's series and edge cases", {
  expect_lt(abs(exp(EDE:::solow1993b_log_fisher(153 / 1714, 28)) -
                  EDE:::solow1993b_fisher_series(153, 1714, 28)), 1e-12)
  expect_equal(exp(EDE:::solow1993b_log_fisher(0.6, 2)), 0.2)
  expect_equal(EDE:::solow1993b_log_fisher(1, 5), 0)
  expect_equal(EDE:::solow1993b_log_fisher(0.2, 5), -Inf)
  # undefined estimate: the other sightings all at the origin
  expect_true(is.na(EDE:::solow1993b_mvue(10, 10, 3)))
})
