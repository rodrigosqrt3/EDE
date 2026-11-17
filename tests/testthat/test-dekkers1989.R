rec <- function(years) sighting_data(data.frame(year = years, sightings = 1))

test_that("moment estimator recovers the index and the end-point", {
  set.seed(21)
  est <- replicate(100, {
    x <- sort(runif(2000), decreasing = TRUE)   # gamma = -1, end-point 1
    f <- EDE:::moment_fit(x, 150)
    c(f$gamma, f$endpoint)
  })
  expect_equal(mean(est[1, ]), -1, tolerance = 0.1)
  expect_equal(median(est[2, ]), 1, tolerance = 0.02)
  # the empirical sd matches the asymptotic one
  expect_equal(sd(est[1, ]), sqrt(EDE:::moment_var(-1) / 150), tolerance = 0.2)
})

test_that("index is about 0 for an exponential tail and positive for Pareto", {
  set.seed(22)
  g0 <- mean(replicate(60, EDE:::moment_fit(sort(1 + rexp(2000), TRUE), 150)$gamma))
  g1 <- mean(replicate(60, EDE:::moment_fit(sort(runif(2000)^-0.5, TRUE), 150)$gamma))
  expect_equal(g0, 0, tolerance = 0.05)
  expect_equal(g1, 0.5, tolerance = 0.05)
})

test_that("asymptotic variance matches the known special values", {
  expect_equal(EDE:::moment_var(0), 1)
  expect_equal(EDE:::moment_var(0.5), 1.25)
  expect_equal(EDE:::moment_var(-1),
               4 * 3 * 8 / (4 * 5))
})

test_that("dekkers1989 returns an ede_estimate with the extra components", {
  set.seed(23)
  sd <- rec(sort(1900 + 100 * (1 - runif(60)^(1 / 3))))
  r <- suppressWarnings(dekkers1989(sd, k = 15))
  expect_s3_class(r, "ede_estimate")
  expect_true(all(c("gamma", "gamma_se", "gamma_ci", "k", "origin") %in%
                    names(r)))
  expect_equal(r$k, 15L)
  expect_equal(r$origin, min(sd$time))
  expect_lt(r$gamma_ci[1], r$gamma)
  expect_gt(r$gamma_ci[2], r$gamma)
})

test_that("warns and returns NA when the index is not negative", {
  set.seed(24)
  sd <- rec(sort(1900 + 10 * rexp(100)))
  expect_warning(r <- dekkers1989(sd, k = 40), "no finite")
  expect_true(is.na(r$estimate))
})

test_that("the origin matters and is validated", {
  sd <- sighting_data(dodo)
  a <- suppressWarnings(dekkers1989(sd, k = 6, origin = 1500))
  b <- suppressWarnings(dekkers1989(sd, k = 6, origin = 1550))
  expect_false(isTRUE(all.equal(a$gamma, b$gamma)))
  expect_error(dekkers1989(sd, k = 6, origin = "a"), "origin")
  # an origin at or after the sightings leaves no usable threshold
  expect_error(dekkers1989(sd, k = 6, origin = 1900), "at least")
  expect_error(dekkers1989(sd, k = 8, origin = min(sd$time)), "cannot exceed")
})

test_that("input validation", {
  sd <- sighting_data(dodo)
  expect_error(dekkers1989(1), "sighting_data")
  expect_error(dekkers1989(sd, k = 2), "at least 3")
  expect_error(dekkers1989(sd, k = 50), "cannot exceed")
  expect_error(dekkers1989(sd, alpha = 2), "alpha")
})

test_that("ede_evi is a data frame over k, with print and plot", {
  sd <- sighting_data(dodo)
  r <- suppressWarnings(ede_evi(sd))
  expect_s3_class(r, "ede_evi")
  expect_named(r, c("k", "gamma", "se", "lower", "upper", "endpoint"))
  expect_equal(r$k, seq(min(r$k), max(r$k)))
  one <- suppressWarnings(dekkers1989(sd, k = r$k[2]))
  expect_equal(r$gamma[2], one$gamma)
  expect_output(print(r), "extreme-value index")
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(plot(r))
})

test_that("a flat top and the internal moment fit are handled", {
  # every latest sighting at the same time (counted as events)
  flat <- sighting_data(data.frame(time = c(0, 5), count = c(0, 10)))
  expect_error(dekkers1989(flat, k = 3, events = TRUE), "span more than one")
  # internal guards of the moment fit
  expect_error(EDE:::moment_fit(c(3, 2, -1), 2), "later than `origin`")
  fit <- EDE:::moment_fit(c(5, 5, 5, 5), 2)
  expect_true(is.na(fit$gamma) && is.na(fit$endpoint))
})
