test_that("jaric2010 reproduces equations 4 and 6", {
  d <- sighting_data(data.frame(
    time = c(1, 3, 6, 10),
    count = 1
  ))

  trend_curve <- jaric2010(d, test_year = 14, data_out = TRUE)
  expect_equal(trend_curve$chance[trend_curve$time == 14], 0.5)
  expect_equal(attr(trend_curve, "average_interval"), 3)
  expect_equal(attr(trend_curve, "trend_coefficient"), 1)

  plain_curve <- jaric2010(d, test_year = 13, data_out = TRUE, trend = FALSE)
  expect_equal(plain_curve$chance[plain_curve$time == 13], 0.5)
})

test_that("jaric2010 returns the upper bound of equations 7 and 8", {
  d <- sighting_data(data.frame(time = c(1, 3, 6, 10), count = 1))
  res <- jaric2010(d, test_year = 100)
  expect_true(is.na(res$estimate))
  expect_equal(res$interval_type, "one-sided")
  expect_equal(res$upper, 10 + 4 * 19)
  expect_equal(res$average_interval, 3)
  expect_equal(res$trend_coefficient, 1)
  expect_equal(jaric2010(d, test_year = 100, trend = FALSE)$upper, 10 + 3 * 19)
})

test_that("jaric2010 reproduces the black-footed ferret example (Table 2)", {
  # months, with June 1972 = 1 and December 1990 = 223
  year <- c(rep(1972, 4), rep(1973, 6), rep(1974, 2), rep(1975, 3),
            rep(1976, 3), 1977, 1979, 1981, 1981, rep(1982, 3), 1983, 1984, 1984)
  month <- c(6, 7, 8, 10, 5:10, 6, 7, 5, 8, 10, 5, 9, 10, 6, 6, 9, 10, 2, 3, 7, 7, 7, 9)
  d <- sighting_data(data.frame(time = (year - 1972) * 12 + month - 5, count = 1))

  plain <- jaric2010(d, test_year = 223, data_out = TRUE, trend = FALSE)
  expect_equal(attr(plain, "average_interval"), 147 / 27)
  expect_lt(abs(plain$chance[plain$time == 223] - 0.068), 1e-3)
  trend <- jaric2010(d, test_year = 223, data_out = TRUE)
  expect_equal(attr(trend, "trend_coefficient"), 1 / 26)
  expect_lt(abs(trend$chance[trend$time == 223] - 0.068), 1e-3)

  # equations 7 and 8 give May and June 1993 (months 252 and 253)
  expect_equal(ceiling(jaric2010(d, test_year = 400, trend = FALSE)$upper), 252)
  expect_equal(ceiling(jaric2010(d, test_year = 400)$upper), 253)
})

test_that("jaric2010 validates trend and minimum records", {
  two <- sighting_data(data.frame(time = c(1, 3), count = 1))
  expect_error(jaric2010(two, alpha = 0, test_year = 20), "must be in")
  expect_error(jaric2010(two), "test_year")
  expect_error(jaric2010(two, test_year = "20"), "test_year")
  expect_error(jaric2010(two, test_year = 20), "at least 3")
  expect_s3_class(jaric2010(two, test_year = 100, trend = FALSE), "ede_estimate")
  expect_error(jaric2010(two, test_year = 20, trend = NA), "TRUE or FALSE")

  nonpositive <- sighting_data(data.frame(time = c(1, 10, 11), count = 1))
  expect_error(jaric2010(nonpositive, test_year = 20), "must be positive")
})

test_that("jaric2010 validates the test horizon and handles no crossing", {
  d <- sighting_data(data.frame(time = c(1, 3, 6, 10), count = 1))
  expect_error(jaric2010(d, test_year = 10), "later than")

  expect_warning(
    res <- jaric2010(d, alpha = 0.01, test_year = 11),
    "never falls"
  )
  expect_true(is.na(res$upper))
})
