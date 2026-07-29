test_that("sighting_data validates and sorts", {
  d <- data.frame(years = c(1920, 1900, 1910), sightings = c(1, 2, 1))
  sd <- sighting_data(d)
  expect_s3_class(sd, "sighting_data")
  expect_equal(sd$time, c(1900, 1910, 1920))
})

test_that("sighting_data rejects bad input", {
  expect_error(sighting_data(data.frame(a = c(1, 1), b = c(1, 2))), "duplicated")
  expect_error(sighting_data(data.frame(a = c(-1, 1), b = c(1, 2))), "negative")
  expect_error(sighting_data(data.frame(a = c("x", "y"), b = c(1, 2))), "numeric")
})

test_that("sighting_data requires `data`", {
  expect_error(sighting_data(), "required")
})

test_that("sighting_data requires a data.frame or matrix", {
  expect_error(sighting_data(list(1, 2)), "data.frame or matrix")
})

test_that("sighting_data requires at least two columns", {
  expect_error(sighting_data(data.frame(a = 1:3)), "at least two columns")
})

test_that("sighting_data rejects NA in time or count", {
  expect_error(
    sighting_data(data.frame(a = c(1, NA), b = c(1, 2))),
    "cannot contain NA"
  )
})

test_that("sighting_data warns when counts exceed all time values", {
  expect_warning(
    sighting_data(data.frame(a = c(1, 2, 3), b = c(1, 2, 100))),
    "check that"
  )
})

test_that("sighting_data accepts a matrix", {
  m <- matrix(c(1900, 1910, 1, 1), ncol = 2)
  sd <- sighting_data(m)
  expect_s3_class(sd, "sighting_data")
  expect_equal(sd$time, c(1900, 1910))
})

test_that("expand_record rejects a record with no positive sightings", {
  sd0 <- sighting_data(data.frame(a = c(1900, 1910), b = c(0, 0)))
  expect_error(solow1993(sd0, 0.05, test_year = 1920), "no sightings")
})

test_that("expand_record rejects end_time before the first sighting", {
  sd1 <- sighting_data(data.frame(a = c(1900, 1910), b = c(1, 1)))
  expect_error(solow1993(sd1, 0.05, test_year = 1890), "cannot precede")
})
