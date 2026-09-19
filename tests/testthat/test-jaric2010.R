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

test_that("jaric2010 returns the first alpha crossing", {
  d <- sighting_data(data.frame(time = c(1, 3, 6, 10), count = 1))
  expect_equal(jaric2010(d, test_year = 100)$estimate, 86)
  expect_equal(jaric2010(d, test_year = 100, trend = FALSE)$estimate, 67)
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
  expect_true(is.na(res$estimate))
})
