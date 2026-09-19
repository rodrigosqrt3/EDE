test_that("print.ede_estimate prints estimate and CI when both are defined", {
  d <- data.frame(years = c(1900, 1905, 1910), sightings = c(1, 1, 1))
  res <- ole(sighting_data(d), alpha = 0.05)
  out <- capture.output(print(res))
  expect_match(out[1], "^<OLE")
  expect_match(paste(out, collapse = " "), "estimate:")
  expect_match(paste(out, collapse = " "), "% CI:")
})

test_that("print.ede_estimate prints Robson estimate and CI", {
  d <- data.frame(
    years = c(1907, 1910, 1915, 1916, 1920, 1925, 1930, 1931),
    sightings = c(1, 1, 3, 4, 3, 1, 2, 1)
  )
  res <- robson1964(sighting_data(d), alpha = 0.05)
  out <- capture.output(print(res))
  expect_length(out, 3L)
  expect_match(out[1], "^<Robson")
  expect_match(out[2], "^  estimate:")
  expect_match(out[3], "one-sided CI:")
})

test_that("print.ede_estimate returns its argument invisibly", {
  d <- data.frame(years = c(1900, 1910), sightings = c(1, 1))
  res <- robson1964(sighting_data(d), alpha = 0.05)
  ret <- withVisible(print(res))
  expect_false(ret$visible)
  expect_identical(ret$value, res)
})

test_that("print.ede_estimate prints a standalone upper confidence bound", {
  d <- sighting_data(data.frame(year = c(1900, 1910, 1920), sightings = 1))
  res <- jaric_roberts2014(d, reliability = c(1, 0.8, 0.5))
  out <- capture.output(print(res))
  expect_match(out[3], "upper confidence bound")
})
