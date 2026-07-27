test_that("print.ede_estimate prints estimate and CI when both are defined", {
  d <- data.frame(years = c(1900, 1905, 1910), sightings = c(1, 1, 1))
  res <- ole(sighting_data(d), alpha = 0.05)
  out <- capture.output(print(res))
  expect_match(out[1], "^<OLE")
  expect_match(paste(out, collapse = " "), "estimate:")
  expect_match(paste(out, collapse = " "), "% CI:")
})

test_that("print.ede_estimate prints estimate only when no CI is defined", {
  d <- data.frame(
    years = c(1907, 1910, 1915, 1916, 1920, 1925, 1930, 1931),
    sightings = c(1, 1, 3, 4, 3, 1, 2, 1)
  )
  res <- robson1964(sighting_data(d), alpha = 0.05)
  out <- capture.output(print(res))
  expect_length(out, 2L)
  expect_match(out[1], "^<Robson")
  expect_match(out[2], "^  estimate:")
})

test_that("print.ede_estimate returns its argument invisibly", {
  d <- data.frame(years = c(1900, 1910), sightings = c(1, 1))
  res <- robson1964(sighting_data(d), alpha = 0.05)
  ret <- withVisible(print(res))
  expect_false(ret$visible)
  expect_identical(ret$value, res)
})
