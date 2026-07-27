test_that("robson1964 reproduces the reference point estimate when the last gap is 1", {
  d <- data.frame(
    years = c(1907, 1910, 1915, 1916, 1920, 1925, 1930, 1931),
    sightings = c(1, 1, 3, 4, 3, 1, 2, 1)
  )
  res <- robson1964(sighting_data(d), alpha = 0.05)
  expect_equal(res$estimate, 1950)
})

test_that("robson1964 uses the real gap between the last two sightings", {
  # this is a regression test for a spacing bug: an earlier draft of this
  # estimator computed the gap between the last two sighting times as a
  # hardcoded 1 instead of the actual spacing, which only happened to be
  # correct when consecutive sightings were exactly 1 time unit apart
  d <- data.frame(years = c(1900, 1905, 1910, 1920, 1935), sightings = c(2, 1, 3, 1, 1))
  res <- robson1964(sighting_data(d), alpha = 0.05)

  expected <- 1935 + 15 * (1 - 0.05) / 0.05
  expect_equal(res$estimate, expected)
  expect_false(isTRUE(all.equal(res$estimate, 1954)))
})

test_that("robson1964 validates alpha and minimum n", {
  d <- data.frame(years = c(1900, 1910), sightings = c(1, 1))
  sd <- sighting_data(d)
  expect_error(robson1964(sd, alpha = 0), "must be in \\(0, 1\\)")
  expect_error(robson1964(sd, alpha = 1), "must be in \\(0, 1\\)")

  d1 <- data.frame(years = 1900, sightings = 1)
  expect_error(robson1964(sighting_data(d1), 0.05), "at least 2")
})
