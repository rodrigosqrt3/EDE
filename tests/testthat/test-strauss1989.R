test_that("strauss1989 reproduces the reference point estimate", {
  d <- data.frame(
    years = c(1907, 1910, 1915, 1916, 1920, 1925, 1930, 1931),
    sightings = c(1, 1, 3, 4, 3, 1, 2, 1)
  )
  res <- strauss1989(sighting_data(d), alpha = 0.05)
  expect_equal(res$estimate, 1943.819, tolerance = 1e-3)
})

test_that("strauss1989 validates alpha and minimum n", {
  d <- data.frame(years = c(1900, 1910), sightings = c(1, 1))
  sd <- sighting_data(d)
  expect_error(strauss1989(sd, alpha = 0), "must be in \\(0, 1\\)")
  expect_error(strauss1989(sd, alpha = 1), "must be in \\(0, 1\\)")

  d1 <- data.frame(years = 1900, sightings = 1)
  expect_error(strauss1989(sighting_data(d1), 0.05), "at least 2")
})

test_that("strauss1989_curve requires at least 2 sightings", {
  d1 <- data.frame(years = 1900, sightings = 1)
  expect_error(strauss1989_curve(sighting_data(d1)), "at least 2")
})

test_that("strauss1989_curve reproduces the point estimate at matching alpha", {
  d <- data.frame(
    years = c(1907, 1910, 1915, 1916, 1920, 1925, 1930, 1931),
    sightings = c(1, 1, 3, 4, 3, 1, 2, 1)
  )
  sd <- sighting_data(d)
  curve <- strauss1989_curve(sd)

  expect_equal(nrow(curve), 100L)
  expect_named(curve, c("time", "chance"))

  row95 <- curve[abs(curve$chance - 0.95) < 1e-9, ]
  pt <- strauss1989(sd, alpha = 0.05)
  expect_equal(row95$time, pt$estimate, tolerance = 1e-3)
})
