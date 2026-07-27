test_that("mcinerny2006 reproduces published values from Table 1 (McInerny et al. 2006)", {
  # each case reconstructs a sighting record with the same n (distinct
  # sighting times), tn, and T as a real row of Table 1 in the paper; the
  # exact placement of intermediate sightings does not affect the formula,
  # only the count and the two endpoints do
  d1 <- sighting_data(data.frame(year = c(2000, 2002, 2003), sightings = c(1, 1, 1)))
  r1 <- mcinerny2006(d1, alpha = 0.05, test_year = 2004, data_out = TRUE)
  expect_equal(r1$chance[r1$time == 2004], 0.333, tolerance = 5e-3)

  d2 <- sighting_data(data.frame(year = c(2000, 2003, 2005), sightings = c(1, 1, 1)))
  r2 <- mcinerny2006(d2, alpha = 0.05, test_year = 2007, data_out = TRUE)
  expect_equal(r2$chance[r2$time == 2007], 0.360, tolerance = 5e-3)

  d3 <- sighting_data(data.frame(
    year = c(2000, 2005, 2010, 2015, 2020, 2025, 2029),
    sightings = rep(1, 7)
  ))
  r3 <- mcinerny2006(d3, alpha = 0.05, test_year = 2030, data_out = TRUE)
  expect_equal(r3$chance[r3$time == 2030], 0.793, tolerance = 5e-3)
})

test_that("mcinerny2006 validates alpha, test_year, and minimum n", {
  d <- sighting_data(data.frame(year = c(1900, 1910), sightings = c(1, 1)))
  expect_error(mcinerny2006(d, alpha = 0, test_year = 2000), "must be in \\(0, 1\\)")
  expect_error(mcinerny2006(d, alpha = 1, test_year = 2000), "must be in \\(0, 1\\)")
  expect_error(mcinerny2006(d, 0.05), "test_year")
  expect_error(mcinerny2006(d, 0.05, test_year = "2000"), "test_year")

  d1 <- sighting_data(data.frame(year = 1900, sightings = 1))
  expect_error(mcinerny2006(d1, 0.05, test_year = 2000), "at least 2")
})

test_that("mcinerny2006 rejects test_year at or before the last sighting", {
  d <- sighting_data(data.frame(year = c(1900, 1910), sightings = c(1, 1)))
  expect_error(mcinerny2006(d, 0.05, test_year = 1910), "later than")
  expect_error(mcinerny2006(d, 0.05, test_year = 1905), "later than")
})

test_that("mcinerny2006 warns and returns NA when persistence is never rejected", {
  # sighting rate is low enough here that chance stays above alpha for a
  # single candidate year right after the last sighting
  d <- sighting_data(data.frame(year = c(1900, 1902, 1904, 1990), sightings = c(1, 1, 1, 1)))
  expect_warning(res <- mcinerny2006(d, alpha = 0.001, test_year = 1991), "never falls")
  expect_true(is.na(res$estimate))
})

test_that("mcinerny2006 warns when the sighting rate is degenerate (n_used / tn >= 1)", {
  # a sighting in every single time unit up to the last one
  d <- sighting_data(data.frame(year = c(1900, 1901, 1902, 1903), sightings = rep(1, 4)))
  expect_warning(mcinerny2006(d, alpha = 0.05, test_year = 1904), "degenerate")
})

test_that("mcinerny2006 point estimate matches the first crossing in its own curve", {
  d <- sighting_data(data.frame(year = c(1900, 1902, 1903, 1905, 1907), sightings = rep(1, 5)))
  curve <- mcinerny2006(d, alpha = 0.05, test_year = 2000, data_out = TRUE)
  pt <- mcinerny2006(d, alpha = 0.05, test_year = 2000)
  expect_equal(pt$estimate, curve$time[curve$chance <= 0.05][1])
})
