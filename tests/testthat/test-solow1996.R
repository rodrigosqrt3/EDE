test_that("two taxa: the statistic is the p-value (Solow 1996)", {
  # Paper: m = 16, X1 = 800, Xm = 1147; n = 8, Y1 = 979, Yn = 1091;
  # A1 = 0.059, and a 0.95 interval for the common end-point of (1147, 1171).
  r <- solow1996(seymour_ammonites[c(3, 8), ])
  expect_equal(round(r$p_value, 3), 0.059)
  expect_equal(exp(-r$statistic / 2), r$p_value, tolerance = 1e-10)
  expect_equal(r$df, 2)
  expect_equal(r$lower, 1147)
  expect_equal(round(r$upper), 1171)
  expect_s3_class(r, "ede_estimate")
  expect_equal(r$interval_type, "one-sided")
})

test_that("three taxa reproduce the paper", {
  # Paper: A2 = 0.013, -2 log A2 = 8.62, p = 0.071, interval (1165, 1182).
  r <- solow1996(seymour_ammonites[c(2, 3, 8), ])
  expect_equal(round(r$statistic, 2), 8.62)
  expect_equal(r$df, 4)
  expect_equal(round(r$p_value, 3), 0.071)
  expect_equal(r$lower, 1165)
  expect_equal(round(r$upper), 1182)
})

test_that("the ten uppermost ammonite taxa reproduce the paper", {
  # Paper: -2 log A = 22.7 on 18 df, p = 0.20; 0.90 and 0.95 intervals
  # (1165, 1174) and (1165, 1178).
  r <- solow1996(seymour_ammonites)
  expect_equal(round(r$statistic, 1), 22.7)
  expect_equal(r$df, 18)
  expect_equal(round(r$p_value, 2), 0.20)
  # the paper rounds its limits to whole meters, not consistently
  expect_lt(abs(solow1996(seymour_ammonites, alpha = 0.1)$upper - 1174), 1)
  expect_lt(abs(r$upper - 1178), 1)
  expect_equal(r$lower, 1165)
  # leaving out taxon 8: the 0.90 upper limit increases to 1183
  r8 <- solow1996(seymour_ammonites[-8, ], alpha = 0.1)
  expect_lt(abs(r8$upper - 1183), 1)
})

test_that("a common lower end-point is tested by changing the sign", {
  # Paper: taxa 3 and 8, common lower end-point: A1 = 0.0013.
  # The common end-point is rejected, so the interval is empty and the
  # function says so.
  expect_warning(
    r <- solow1996(seymour_ammonites[c(3, 8), ], end = "lower"), "empty"
  )
  expect_equal(round(r$p_value, 4), 0.0013)
  expect_lte(r$lower, r$upper)
  expect_equal(r$upper, 800)
})

test_that("one taxon gives the interval of Strauss & Sadler", {
  d <- seymour_ammonites[3, ]
  r <- solow1996(d)
  expect_true(is.na(r$p_value))
  expect_equal(r$df, 0)
  sd <- sighting_data(data.frame(time = seq(800, 1147, length.out = 16),
                                 count = 1))
  s <- strauss1989(sd)
  expect_equal(r$upper, s$upper, tolerance = 1e-7)
  expect_output(print(r), "not defined for one taxon")
})

test_that("a list of sighting_data objects is equivalent to the summary", {
  a <- sighting_data(data.frame(time = c(1, 3, 4, 6, 9), count = c(1, 2, 1, 1, 1)))
  b <- sighting_data(data.frame(time = c(2, 5, 8), count = 1))
  r1 <- solow1996(list(a, b))
  r2 <- solow1996(data.frame(n = c(5, 3), lowest = c(1, 2), highest = c(9, 8)))
  expect_equal(r1$upper, r2$upper)
  expect_equal(r1$p_value, r2$p_value)
  expect_output(print(r1), "common upper end-point")
})

test_that("input is validated", {
  expect_error(solow1996(1), "data frame")
  expect_error(solow1996(sighting_data(data.frame(t = 1:3, c = 1))), "one entry")
  expect_error(solow1996(data.frame(n = 1, lowest = 0, highest = 1)),
               "at least 2")
  expect_error(solow1996(data.frame(n = 5, lowest = 2, highest = 1)),
               "greater than")
  expect_error(solow1996(data.frame(n = 5, low = 2, high = 1)), "columns")
  expect_error(solow1996(seymour_ammonites, alpha = 0), "alpha")
  # a record rejected at the highest find returns it, with a warning
  strong <- data.frame(n = c(50, 50), lowest = c(0, 0), highest = c(100, 30))
  expect_warning(r <- solow1996(strong), "empty")
  expect_equal(r$upper, 100)
})

test_that("a taxon without sightings and non-finite input are refused", {
  ok <- sighting_data(data.frame(time = 1:3, count = 1))
  none <- sighting_data(data.frame(time = 1:3, count = 0))
  expect_error(solow1996(list(ok, none)), "no sightings")
  expect_error(solow1996(data.frame(n = 5, lowest = NA_real_, highest = 3)),
               "finite numbers")
})
