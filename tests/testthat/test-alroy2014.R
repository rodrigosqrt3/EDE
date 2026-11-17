# The illustrative example of the paper: a range of r = 50 time units with
# ps = 0.25: 14 presences (first and last included) and 36 absences, so that
# (np - 2) / (np - 2 + na) = 12 / 48.
toy_times <- c(0, 4, 7, 11, 15, 18, 22, 26, 29, 33, 37, 41, 45, 49)
toy <- sighting_data(data.frame(time = toy_times, count = 1))

test_that("the worked example of the paper is reproduced", {
  # Paper: lambda = 0.0139, E = 0.0138, e1 = 0.0183, e = 0.072 after three
  # intervals, 0.506 after ten, and past 0.99 after 25.
  s <- alroy2014(toy, test_year = 49 + 40, adjust_range = FALSE,
                 data_out = TRUE)
  expect_equal(round(s$e[1], 4), 0.0183)
  expect_equal(round(s$e[3], 3), 0.072)
  expect_equal(round(s$e[10], 3), 0.506)
  expect_gt(s$e[25], 0.989)
  expect_equal(which(s$e > 0.99)[1], 25)
  expect_equal(s$time[1], 50)
  res <- alroy2014(toy, test_year = 89, adjust_range = FALSE)
  expect_equal(res$r, 50)
  expect_equal(res$n_presence, 14)
  expect_equal(res$n_absence, 36)
  expect_equal(res$ps, 0.25)
  expect_equal(-log(0.5) / 50, -log(1 - res$E), tolerance = 1e-10)
  expect_equal(round(res$E, 4), 0.0138)
})

test_that("the recursion is equation 4 and collapses to equation 3", {
  res <- alroy2014(toy, test_year = 200, adjust_range = FALSE)
  e <- res$shadow$e
  big_e <- res$E
  ps <- res$ps
  expect_equal(e[1], big_e / (big_e + (1 - big_e) * (1 - ps)))
  a <- e[1] + (1 - e[1]) * big_e
  expect_equal(e[2], a / (a + (1 - e[1]) * (1 - big_e) * (1 - ps)))
  expect_true(all(diff(e[1:60]) > 0)) # increasing until it saturates at 1
  expect_true(all(e > 0 & e <= 1))
})

test_that("the range is corrected for sighting failures by default", {
  adj <- alroy2014(toy, test_year = 300)
  raw <- alroy2014(toy, test_year = 300, adjust_range = FALSE)
  r_adj <- 50 * (14 + 1) / (14 - 1)
  expect_equal(adj$E, 1 - 0.5^(1 / r_adj))
  expect_lt(adj$E, raw$E) # a longer range means a lower extinction rate
  expect_gt(adj$shadow$e[10], 0)
  expect_lt(adj$shadow$e[10], raw$shadow$e[10])
})

test_that("the run-based probabilities reproduce the example of the paper", {
  # Paper: p-a-a-p-p-a-a-p gives ps|p = 1/3 and ps|a = 2/4.
  sd <- sighting_data(data.frame(time = c(0, 3, 4, 7), count = 1))
  res <- alroy2014(sd, test_year = 200, run_based = TRUE)
  expect_equal(unname(res$ps["first"]), 1 / 3)
  expect_equal(unname(res$ps["next_unit"]), 2 / 4)
  # the first shadow value uses ps|p and later ones ps|a
  s <- res$shadow$e
  big_e <- res$E
  expect_equal(s[1], big_e / (big_e + (1 - big_e) * (1 - 1 / 3)))
  a <- s[1] + (1 - s[1]) * big_e
  expect_equal(s[2], a / (a + (1 - s[1]) * (1 - big_e) * (1 - 0.5)))
  # without embedded absences or presences the paper's convention applies
  none <- sighting_data(data.frame(time = c(0, 1, 2, 3), count = 1))
  r0 <- suppressWarnings(alroy2014(none, test_year = 6, run_based = TRUE))
  expect_equal(unname(r0$ps["first"]), 3 / (3 + 1)) # npp / (npp + 1), npa = 0
  expect_equal(unname(r0$ps["next_unit"]), 0 / (0 + 1)) # nap / (nap + 1)
})

test_that("a record with no absences gives certainty at once", {
  full <- sighting_data(data.frame(time = 0:5, count = 1))
  s <- alroy2014(full, test_year = 8, data_out = TRUE)
  expect_equal(s$e, c(1, 1, 1))
})

test_that("estimate and credible limits follow the shadow values", {
  s <- alroy2014(toy, test_year = 90, adjust_range = FALSE, data_out = TRUE)
  one <- alroy2014(toy, test_year = 90, adjust_range = FALSE)
  expect_equal(one$estimate, s$time[which(s$e >= 0.5)[1]]) # 50 + 9 = 59
  expect_equal(one$estimate, 59)
  expect_equal(one$upper, s$time[which(s$e >= 0.95)[1]])
  expect_true(is.na(one$lower))
  expect_equal(one$interval_kind, "credible")
  expect_equal(one$interval_type, "one-sided")
  two <- alroy2014(toy, test_year = 90, adjust_range = FALSE,
                   interval = "two-sided", alpha = 0.5)
  expect_equal(two$lower, s$time[which(s$e >= 0.25)[1]])
  expect_equal(two$upper, s$time[which(s$e >= 0.75)[1]])
  expect_s3_class(one, "ede_estimate")
  expect_output(print(one), "Alroy")
})

test_that("a limit that is not reached gives NA with a warning", {
  expect_warning(
    res <- alroy2014(toy, test_year = 52, adjust_range = FALSE),
    "does not reach"
  )
  expect_true(is.na(res$estimate))
  expect_true(is.na(res$upper))
})

test_that("a truncated range and an earlier base change the result", {
  base_res <- alroy2014(toy, test_year = 400, adjust_range = FALSE)
  trunc <- alroy2014(toy, test_year = 400, adjust_range = FALSE,
                     truncated = TRUE)
  expect_equal(trunc$ps, (14 - 1) / (14 - 1 + 36))
  expect_gt(trunc$ps, base_res$ps)
  early <- alroy2014(toy, test_year = 400, adjust_range = FALSE, base = -50)
  expect_equal(early$r, 100)
  # a longer range gives a lower rate and a later extinction estimate
  expect_lt(early$E, base_res$E)
  expect_gt(early$estimate, base_res$estimate)
})

test_that("input is validated", {
  expect_error(alroy2014(toy, test_year = 40), "later than the last")
  expect_error(alroy2014(toy), "test_year")
  expect_error(alroy2014(toy, test_year = 60.5), "whole number")
  expect_error(alroy2014(toy, test_year = 60, base = 10), "not later")
  expect_error(alroy2014(toy, test_year = 60, alpha = 1), "alpha")
  few <- sighting_data(data.frame(time = c(0, 5), count = 1))
  expect_error(alroy2014(few, test_year = 10), "at least 3")
  frac <- sighting_data(data.frame(time = c(0, 1.5, 4), count = 1))
  expect_error(alroy2014(frac, test_year = 10), "whole numbers")
})
