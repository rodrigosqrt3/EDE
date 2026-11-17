rec <- function(x) sighting_data(data.frame(year = x, sightings = 1))

test_that("profile log-likelihood is maximised by the reported estimate", {
  set.seed(31)
  all <- sort(1 - runif(300)^(1 / 4), decreasing = TRUE)
  x <- all[1:40]; u <- all[41]
  fit <- EDE:::smith_fit(x, u, x[1] - u)
  expect_false(is.na(fit$s))
  g <- x[1] - u
  nb <- vapply(fit$s * c(0.9, 1.1), EDE:::smith_profile, numeric(1),
               x = x, u = u, g = g)
  expect_true(all(nb < fit$ll))
})

test_that("closed-form shape equals the numerical maximiser over the shape", {
  set.seed(32)
  all <- sort(1 - runif(200)^(1 / 3), decreasing = TRUE)
  x <- all[1:30]; u <- all[31]; s <- 0.05
  r <- (x[1] - x + s) / (x[1] - u + s)
  ll_a <- function(a) 30 * log(a) - 30 * log(x[1] - u + s) + (a - 1) * sum(log(r))
  a_hat <- -30 / sum(log(r))
  expect_equal(optimize(ll_a, c(0.01, 100), maximum = TRUE)$maximum, a_hat,
               tolerance = 1e-4)
  expect_equal(EDE:::smith_profile(s, x, u, x[1] - u), ll_a(a_hat),
               tolerance = 1e-8)
})

test_that("smith1985 estimates the end-point in the regular range", {
  set.seed(33)
  r <- replicate(40, {
    sd <- rec(sort(1 - runif(400)^(1 / 4)))
    out <- suppressWarnings(smith1985(sd, k = 40))
    c(out$estimate, out$shape)
  })
  ok <- !is.na(r[1, ])
  expect_gt(mean(ok), 0.7)
  expect_equal(median(r[1, ok]), 1, tolerance = 0.2)
  expect_gt(median(r[2, ok]), 2)
})

test_that("warns when the problem is not regular", {
  expect_warning(
    expect_warning(smith1985(sighting_data(dodo), k = 8), "not larger than 2"),
    "upper confidence limit")
})

test_that("returns NA with a warning when there is no interior maximum", {
  sd <- rec(c(1900, 1904, 1909, 1911, 1913, 1918, 1920, 1925, 1928, 1930, 1934))
  expect_warning(r <- smith1985(sd, k = 6), "no interior maximum")
  expect_true(is.na(r$estimate))
})

test_that("components and validation", {
  r <- suppressWarnings(smith1985(sighting_data(dodo), k = 8))
  expect_s3_class(r, "ede_estimate")
  expect_equal(r$k, 8L)
  expect_equal(r$gamma, -1 / r$shape)
  expect_false(r$regular)
  expect_gte(r$lower, max(sighting_data(dodo)$time[sighting_data(dodo)$count > 0]) -
               1e-6)
  expect_error(smith1985(1), "sighting_data")
  expect_error(smith1985(sighting_data(dodo), k = 3), "at least 4")
  expect_error(smith1985(sighting_data(dodo), alpha = 0), "alpha")
})

test_that("a flat top is refused", {
  flat <- sighting_data(data.frame(time = c(0, 1), count = c(1, 30)))
  expect_error(smith1985(flat, k = 4, events = TRUE), "span more than one")
})

test_that("an interior maximum below the exponential limit is not accepted", {
  testthat::local_mocked_bindings(
    smith_profile = function(s, x, u, g) -(log(s / g))^2,
    smith_profile_inf = function(x, u) 5
  )
  fit <- EDE:::smith_fit(c(10, 9, 8, 7), 6, 4)
  expect_true(is.na(fit$s) && is.na(fit$shape) && is.na(fit$ll))
})

test_that("an upper limit that is never crossed falls back to the grid end", {
  testthat::local_mocked_bindings(
    smith_profile = function(s, x, u, g) 1,
    smith_profile_inf = function(x, u) -1
  )
  ci <- EDE:::smith_interval(c(10, 9, 8, 7), 6, 4,
                             list(s = 2, ll = 1, shape = 3), crit = 1)
  expect_true(is.finite(ci[2]))
  expect_equal(ci[2], 2 * exp(log(1e5 * 4 / 2)))
})
