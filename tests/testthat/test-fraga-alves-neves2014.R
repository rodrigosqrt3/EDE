years <- c(1900, 1904, 1909, 1911, 1913, 1918, 1920, 1925, 1928, 1930, 1934)
rec <- function() sighting_data(data.frame(year = years, sightings = 1))

test_that("fraga_alves_neves2014 matches the closed form of equation 3", {
  k <- 5
  t <- sort(years, decreasing = TRUE)
  a <- (log(k + 0:(k - 1) + 1) - log(k + 0:(k - 1))) / log(2)
  expected <- t[1] + sum(a * (t[k + 1] - t[k + 1 + 0:(k - 1)]))
  res <- fraga_alves_neves2014(rec(), k = k)
  expect_s3_class(res, "ede_estimate")
  expect_equal(res$estimate, expected)
  expect_equal(res$k, 5L)
  # weights are positive and sum to one
  expect_equal(sum(a), 1)
  expect_true(all(a > 0))
})

test_that("fraga_alves_neves2014 equals the form of equation 16", {
  k <- 4
  t <- sort(years, decreasing = TRUE)
  eq16 <- t[k + 1] + t[1] +
    sum(log((k + 0:(k - 1)) / (k + 0:(k - 1) + 1)) * t[k + 1 + 0:(k - 1)]) /
      log(2)
  expect_equal(fraga_alves_neves2014(rec(), k = k)$estimate, eq16)
})

test_that("fraga_alves_neves2014 never falls below the last sighting", {
  res <- fraga_alves_neves2014(rec())
  expect_gt(res$estimate, 1934)
  expect_equal(res$k, 5L)            # floor(11 / 2)
  expect_true(is.na(res$lower) && is.na(res$upper))
  out <- capture.output(print(res))
  expect_length(out, 2L)
  expect_match(out[2], "estimate:")
})

test_that("fraga_alves_neves2014 improves on the sample maximum", {
  skip_on_cran()
  # tail 1 - F(x) = exp(-1 / (1 - x)) lies in the Gumbel domain, end-point 1
  set.seed(10)
  est <- replicate(100, {
    x <- 1 + 1 / log(stats::runif(3000))
    d <- sighting_data(data.frame(time = x, count = 1))
    c(fraga_alves_neves2014(d, k = 150)$estimate, max(x))
  })
  expect_true(all(est[1, ] > est[2, ]))
  expect_lt(sqrt(mean((est[1, ] - 1)^2)), sqrt(mean((est[2, ] - 1)^2)))
})

test_that("fraga_alves_neves2014 counts events only when asked", {
  d <- sighting_data(data.frame(year = c(1900, 1910, 1920, 1930, 1940),
                                sightings = c(1, 2, 1, 2, 1)))
  distinct <- fraga_alves_neves2014(d, k = 2)
  events <- fraga_alves_neves2014(d, k = 3, events = TRUE)
  expect_equal(distinct$k, 2L)
  expect_equal(events$k, 3L)
  expect_error(fraga_alves_neves2014(d, k = 3), "cannot exceed")
})

test_that("fraga_alves_neves2014 validates its inputs", {
  expect_error(fraga_alves_neves2014(data.frame(1)), "sighting_data")
  expect_error(fraga_alves_neves2014(rec(), k = 0), "at least 1")
  expect_error(fraga_alves_neves2014(rec(), k = 6), "cannot exceed 5")
  one <- sighting_data(data.frame(year = 1900, sightings = 1))
  expect_error(fraga_alves_neves2014(one), "at least 2")
  none <- sighting_data(data.frame(year = c(1900, 1910), sightings = 0))
  expect_error(fraga_alves_neves2014(none), "no sightings")
})
