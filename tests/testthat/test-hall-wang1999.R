# Men's 100 m records at the 1988 and 1992 Olympics (Hall & Wang 1999,
# Table 3). The paper estimates a lower end-point; EDE estimates the end-point
# at the largest time, so the record is entered as time = -seconds.
sprint <- function() {
  secs <- c(9.92, 9.96, 9.99, 10.02, 10.04, 10.09, 10.10, 10.11, 10.12, 10.22,
            10.26)
  sighting_data(data.frame(time = -secs, count = ifelse(secs == 10.11, 2, 1)))
}

test_that("simulated Greenwood null matches the exact points of Burrows (1979)", {
  # 10 G for 10 uniform spacings (m = 11 exponentials), quoted by Stephens
  # (1981): lower 1% and 5%, upper 5% and 1%
  q <- stats::quantile(EDE:::greenwood_null(11L, 1e5) * 10,
                       c(0.01, 0.05, 0.95, 0.99), names = FALSE)
  expect_equal(q, c(1.117, 1.211, 2.404, 3.008), tolerance = 0.01)
})

test_that("hall_wang1999 reproduces the 100 m example of Hall & Wang (1999)", {
  res <- hall_wang1999(sprint(), alpha = 0.1, k = 12, events = TRUE)
  expect_s3_class(res, "ede_estimate")
  # paper: minimum of G at 9.91, extremal exponent 0.82, interval (9.80, 9.92)
  expect_equal(-res$estimate, 9.91, tolerance = 0.01)
  expect_equal(res$shape, 0.82, tolerance = 0.01)
  expect_equal(-res$upper, 9.80, tolerance = 0.005)
  expect_equal(-res$lower, 9.92, tolerance = 0.005)
  expect_equal(res$k, 12L)
  expect_equal(res$method, "Hall & Wang (1999)")
})

test_that("hall_wang1999 is location- and scale-equivariant", {
  base <- hall_wang1999(sprint(), alpha = 0.1, k = 12, events = TRUE)
  d <- sprint()
  d$time <- 3 * d$time + 100
  moved <- hall_wang1999(d, alpha = 0.1, k = 12, events = TRUE)
  expect_equal(moved$estimate, 3 * base$estimate + 100, tolerance = 1e-6)
  expect_equal(moved$upper, 3 * base$upper + 100, tolerance = 1e-6)
  expect_equal(moved$shape, base$shape, tolerance = 1e-6)
})

test_that("hall_wang1999 puts the estimate after the last sighting", {
  d <- sighting_data(data.frame(
    year = c(1900, 1904, 1909, 1911, 1913, 1918, 1920, 1925, 1928, 1930, 1934),
    sightings = 1
  ))
  # the default k = 10 gives no upper limit at the 5% level
  res <- suppressWarnings(hall_wang1999(d))
  expect_gt(res$estimate, 1934)
  expect_equal(res$lower, 1934, tolerance = 1e-3)
  expect_true(res$k == 10L)
})

test_that("hall_wang1999 reports a missing minimum as NA with warnings", {
  d <- sighting_data(data.frame(time = c(43, 28, 27, 26, 25), count = 1))
  expect_warning(
    expect_warning(res <- hall_wang1999(d, k = 5), "does not exist for this"),
    "upper confidence limit does not exist"
  )
  expect_true(is.na(res$estimate))
  expect_true(is.infinite(res$upper))
  expect_gt(res$lower, 43)
})

test_that("the p-value curve is consistent with the interval", {
  res <- hall_wang1999(sprint(), alpha = 0.1, k = 12, events = TRUE)
  cv <- hall_wang1999(sprint(), alpha = 0.1, k = 12, events = TRUE,
                      data_out = TRUE, test_year = -9.5)
  expect_named(cv, c("time", "statistic", "p_value"))
  expect_true(all(cv$p_value >= 0 & cv$p_value <= 1))
  inside <- cv$time[cv$p_value > 0.1]
  expect_equal(min(inside), res$lower, tolerance = 0.01)
  expect_equal(max(inside), res$upper, tolerance = 0.01)
  expect_error(hall_wang1999(sprint(), k = 12, events = TRUE, data_out = TRUE),
               "test_year")
})

test_that("the critical values do not disturb the random number stream", {
  set.seed(7)
  expected <- runif(1)
  set.seed(7)
  suppressWarnings(hall_wang1999(sprint(), k = 8, events = TRUE, nsim = 2000))
  expect_equal(runif(1), expected)
  rm(".Random.seed", envir = globalenv())
  suppressWarnings(hall_wang1999(sprint(), k = 9, events = TRUE, nsim = 2000))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("the null distribution of Greenwood's statistic has the exact mean", {
  # E[sum D_i^2] = 2 / (m + 1) for m uniform spacings (Dirichlet(1, ..., 1))
  for (m in c(3L, 11L)) {
    expect_equal(mean(greenwood_null(m, 1e5)), 2 / (m + 1), tolerance = 0.01)
  }
})

test_that("confidence intervals attain their nominal coverage", {
  skip_on_cran()
  set.seed(2024)
  covered <- replicate(300, {
    # Weibull extreme-value model with extremal exponent 2, end-point 0
    gam <- cumsum(stats::rexp(8))^(1 / 2)
    d <- sighting_data(data.frame(time = -gam * 10, count = 1))
    r <- suppressWarnings(hall_wang1999(d, alpha = 0.1, k = 8, nsim = 20000))
    isTRUE(r$lower <= 0 && 0 <= r$upper)
  })
  expect_gt(mean(covered), 0.85)
  expect_lt(mean(covered), 0.96)
})

test_that("hall_wang1999 validates its inputs", {
  d <- sprint()
  expect_error(hall_wang1999(data.frame(1)), "sighting_data")
  expect_error(hall_wang1999(d, alpha = 0), "must be in \\(0, 1\\)")
  expect_error(hall_wang1999(d, alpha = 1), "must be in \\(0, 1\\)")
  expect_error(hall_wang1999(d, k = 3), "at least 4")
  expect_error(hall_wang1999(d, k = 5.5), "at least 4")
  expect_error(hall_wang1999(d, k = 50), "cannot exceed")
  expect_error(hall_wang1999(d, nsim = 10), "at least 1000")
  few <- sighting_data(data.frame(year = c(1900, 1910, 1920), sightings = 1))
  expect_error(hall_wang1999(few), "at least 4")
  none <- sighting_data(data.frame(year = c(1900, 1910), sightings = 0))
  expect_error(hall_wang1999(none), "no sightings")
})

test_that("a flat top is refused and an empty interval is reported", {
  # 30 events at one time: the 4 latest sightings do not span more than one time
  flat <- sighting_data(data.frame(time = c(0, 1), count = c(1, 30)))
  expect_error(hall_wang1999(flat, k = 4, events = TRUE, nsim = 2000),
               "span more than one")
  t <- c(1906, 1916, 1918, 1919, 1924, 1933, 1952, 1957, 1959, 1981, 1982, 1985)
  expect_warning(
    res <- hall_wang1999(sighting_data(data.frame(time = t, count = 1)),
                         alpha = 0.999, k = 6, nsim = 2000),
    "interval is empty"
  )
  expect_true(is.na(res$lower) && is.na(res$upper))
})

test_that("an interval that only starts beyond the grid is returned, not an error", {
  # gaps whose statistic stays above its limit at infinity over the whole
  # grid; a critical value between the two leaves only the limit inside
  g <- c(0, 9.74, 76.34, 96.08, 112.85, 137, 137.63)
  grid <- EDE:::hw_grid(g)
  st <- EDE:::hw_stat(grid, g)
  inf <- EDE:::hw_stat_inf(g)
  expect_true(all(st > inf))
  crit <- (min(st) + inf) / 2
  expect_equal(EDE:::hw_interval(g, NULL, crit), c(grid[length(grid)], Inf))
})
