test_that("solow1993b reproduces published value for Black-footed ferret (Solow 1993b)", {
  p <- EDE:::solow1993b_fisher_series(153, 1714, 28) / EDE:::solow1993b_fisher_series(229, 1714, 28)
  expect_equal(round(p, 3), 0.050)
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
})

test_that("solow1993b estimates extinction year and handles warning when alpha is not reached", {
  d <- sighting_data(data.frame(
    years = c(1907, 1910, 1915, 1916, 1920, 1925, 1930, 1931),
    sightings = c(1, 1, 3, 4, 3, 1, 2, 1)
  ))
  res <- solow1993b(d, alpha = 0.05, test_year = 2000)
  expect_false(is.na(res$estimate))

  expect_warning(
    res_na <- solow1993b(d, alpha = 1e-10, test_year = 1932),
    "never falls"
  )
  expect_true(is.na(res_na$estimate))
})
