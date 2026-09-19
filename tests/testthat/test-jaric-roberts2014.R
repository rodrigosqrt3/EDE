test_that("jaric_roberts2014 implements the reliability-adjusted equations", {
  d <- sighting_data(data.frame(
    year = c(1900, 1910, 1920),
    sightings = 1
  ))
  res <- jaric_roberts2014(d, reliability = c(1, 0.8, 0.5))

  expect_equal(res$effective_sightings, 1.3)
  expect_equal(res$effective_endpoint, 1914)
  expect_equal(res$estimate, 1924.7692308, tolerance = 1e-7)
  expect_equal(res$upper, 2040.2556419, tolerance = 1e-7)
  expect_equal(res$interval_type, "one-sided")
})

test_that("jaric_roberts2014 reduces to the Solow estimates at reliability one", {
  d <- sighting_data(data.frame(
    year = c(1900, 1910, 1920, 1930),
    sightings = 1
  ))
  res <- jaric_roberts2014(d, reliability = rep(1, 4))

  expect_equal(res$effective_sightings, 3)
  expect_equal(res$effective_endpoint, 1930)
  expect_equal(res$estimate, 1900 + (4 / 3) * 30)
  expect_equal(res$upper, 1900 + 30 / 0.05^(1 / 3))
})

test_that("jaric_roberts2014 accepts row-aligned reliability", {
  d <- sighting_data(data.frame(
    year = c(1900, 1905, 1910, 1920),
    sightings = c(0, 1, 0, 1)
  ))
  res <- jaric_roberts2014(d, reliability = c(0, 0.8, 0, 0.5))
  expect_equal(res$effective_sightings, 1.3)
})

test_that("jaric_roberts2014 accepts reliability for occupied rows only", {
  d <- sighting_data(data.frame(
    year = c(1900, 1905, 1910, 1920),
    sightings = c(0, 1, 0, 1)
  ))
  res <- jaric_roberts2014(d, reliability = c(0.8, 0.5))
  expect_equal(res$effective_sightings, 1.3)
})

test_that("jaric_roberts2014 validates reliability and binary counts", {
  d <- sighting_data(data.frame(year = c(1900, 1910), sightings = 1))
  expect_error(jaric_roberts2014(d, c(1, 1), alpha = 0), "must be in")
  expect_error(jaric_roberts2014(d, c(1, 1.2)), "in \\[0, 1\\]")
  expect_error(jaric_roberts2014(d, 1), "one value per")

  counted <- sighting_data(data.frame(year = c(1900, 1910), sightings = c(1, 2)))
  expect_error(jaric_roberts2014(counted, c(1, 1)), "binary")
})

test_that("jaric_roberts2014 requires usable reliability after the origin", {
  origin_only <- sighting_data(data.frame(year = 1900, sightings = 1))
  expect_error(
    jaric_roberts2014(origin_only, reliability = 1),
    "at least one sighting"
  )

  zero_reliability <- sighting_data(data.frame(
    year = c(1900, 1910),
    sightings = 1
  ))
  expect_error(
    jaric_roberts2014(zero_reliability, reliability = c(1, 0)),
    "sum of reliabilities"
  )
})
