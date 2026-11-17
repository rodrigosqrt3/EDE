rec <- function(years) sighting_data(data.frame(year = years, sightings = 1))
sd <- rec(c(1900, 1904, 1909, 1911, 1913, 1918, 1920, 1925))

test_that("equation 2 of McCarthy (1998) for a hand-computed case", {
  effort <- ifelse(1900:1950 < 1915, 1, 2)
  res <- mccarthy1998(sd, effort, test_year = 1950, data_out = TRUE)
  cum <- cumsum(effort)
  n <- 7   # 8 sighting years, the one at the origin is omitted
  e_last <- cum[1925 - 1900 + 1]
  expect_equal(res$p_value[res$time == 1930],
               (e_last / cum[1930 - 1900 + 1])^n)
  expect_equal(res$time, 1926:1950)
  expect_true(all(diff(res$p_value) < 0))
})

test_that("constant effort gives the Solow test up to the discretisation", {
  res <- mccarthy1998(sd, rep(1, 51), test_year = 1950, data_out = TRUE)
  expect_equal(res$p_value[res$time == 1940],
               ((1925 - 1900 + 1) / (1940 - 1900 + 1))^7)
  # the upper bounds are close
  a <- mccarthy1998(sd, rep(1, 101), alpha = 0.05, test_year = 2000)
  b <- solow1993(sd, alpha = 0.05, test_year = 2000)
  expect_lt(abs(a$upper - b$upper), 3)
})

test_that("rising effort makes the bound earlier, falling effort later", {
  base <- mccarthy1998(sd, rep(1, 101), test_year = 2000)$upper
  rising <- mccarthy1998(sd, c(rep(1, 26), rep(5, 75)), test_year = 2000)$upper
  falling <- suppressWarnings(
    mccarthy1998(sd, c(rep(1, 26), rep(0.2, 75)), test_year = 2000)$upper)
  expect_lt(rising, base)
  expect_true(is.na(falling) || falling > base)
})

test_that("effort as a function equals effort as a vector", {
  f <- function(t) 1 + (t - 1900) / 50
  a <- mccarthy1998(sd, f, test_year = 1960)
  b <- mccarthy1998(sd, f(1900:1960), test_year = 1960)
  expect_equal(a$upper, b$upper)
})

test_that("bound is the first time unit at which p <= alpha", {
  eff <- rep(1, 101)
  res <- mccarthy1998(sd, eff, alpha = 0.1, test_year = 2000)
  curve <- mccarthy1998(sd, eff, test_year = 2000, data_out = TRUE)
  expect_equal(res$upper, min(curve$time[curve$p_value <= 0.1]))
  expect_equal(res$alpha, 0.1)
  expect_s3_class(res, "ede_estimate")
  expect_true(is.na(res$estimate))
})

test_that("warns and returns NA when the bound is beyond test_year", {
  expect_warning(r <- mccarthy1998(sd, rep(1, 31), test_year = 1930),
                 "never falls")
  expect_true(is.na(r$upper))
})

test_that("input validation", {
  expect_error(mccarthy1998(1, 1, test_year = 2000), "sighting_data")
  expect_error(mccarthy1998(sd, rep(1, 10), test_year = 1950), "must cover")
  expect_error(mccarthy1998(sd, rep(-1, 51), test_year = 1950),
               "non-negative")
  expect_error(mccarthy1998(sd, rep(NA_real_, 51), test_year = 1950),
               "non-negative")
  expect_error(mccarthy1998(sd, rep(1, 51), test_year = 1920), "later than")
  expect_error(mccarthy1998(sd, rep(1, 51), test_year = 1950.5), "whole")
  expect_error(mccarthy1998(sd, rep(1, 51), alpha = 2, test_year = 1950),
               "alpha")
  expect_error(mccarthy1998(rec(c(1900.5, 1910)), rep(1, 51),
                            test_year = 1950), "whole numbers")
  expect_error(mccarthy1998(sd, c(0, 0, rep(1, 49)), test_year = 1950),
               NA)
  expect_error(mccarthy1998(sd, rep(0, 51), test_year = 1950),
               "no sampling effort")
})

test_that("McCarthy needs a sighting after the observation origin", {
  origin_only <- sighting_data(data.frame(year = 1900, count = 1))
  expect_error(mccarthy1998(origin_only, rep(1, 20), test_year = 1910),
               "after the observation origin")
})
