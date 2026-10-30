# Marshall (2010), "Murder in the Dead Bones": six Geiger-counter clicks at
# 9:15, 9:30, 9:40, 10:05, 10:25 and 10:55 p.m., coded in minutes after 9:15.
dead_bones <- function() {
  sighting_data(data.frame(time = c(0, 15, 25, 50, 70, 100), count = 1))
}
radon224 <- function(t) exp(-0.693 * t / 108)

test_that("constant recovery reproduces the classical intervals (Marshall 2010, section 3)", {
  sd <- dead_bones()
  # 50% interval 14.9 min, 95% interval 82 min after the last click (p. 294)
  expect_equal(round(marshall1997(sd, alpha = 0.5)$upper - 100, 1), 14.9)
  expect_equal(round(marshall1997(sd, alpha = 0.05)$upper - 100), 82)
  # equals strauss1989() for constant recovery
  expect_lt(abs(marshall1997(sd)$upper - strauss1989(sd)$upper), 1e-6)
  # confidence of departure: 5% by 10:56, 73% by 11:25 (equation 3),
  # 38% by 11:05 (section 7)
  expect_equal(round(marshall1997(sd, test_year = 101)$confidence, 2), 0.05)
  expect_equal(round(marshall1997(sd, test_year = 130)$confidence, 2), 0.73)
  expect_equal(round(marshall1997(sd, test_year = 110)$confidence, 2), 0.38)
})

test_that("exponential recovery reproduces Case 2 (radon-224, Marshall 2010, p. 296-298)", {
  sd <- dead_bones()
  # 50% interval ~22 min (to 11:17 p.m.), 95% interval 209 min (to 2:24 a.m.)
  expect_equal(round(marshall1997(sd, radon224, alpha = 0.5)$upper - 100), 22)
  expect_equal(round(marshall1997(sd, radon224, alpha = 0.05)$upper - 100), 209)
  # confidence of departure before 11:25 p.m. drops from 73% to 59%
  expect_equal(round(marshall1997(sd, radon224, test_year = 130)$confidence, 2), 0.59)
})

test_that("a recovery potential that dies out gives an unbounded interval (Case 3)", {
  bismuth211 <- function(t) exp(-0.693 * t / 2.14)
  res <- marshall1997(dead_bones(), bismuth211, alpha = 0.5)
  expect_equal(res$upper, Inf)
})

test_that("the interval is bounded by max_time and scale-invariant in recovery", {
  sd <- dead_bones()
  a <- marshall1997(sd, radon224)
  b <- marshall1997(sd, function(t) 7 * radon224(t))
  expect_lt(abs(a$upper - b$upper), 1e-6)
  expect_equal(marshall1997(sd, max_time = 150)$upper, Inf)
  expect_lt(abs(marshall1997(sd, max_time = 1000)$upper - strauss1989(sd)$upper), 1e-6)
})

test_that("marshall1997 returns a one-sided ede_estimate", {
  res <- marshall1997(dead_bones(), radon224)
  expect_s3_class(res, "ede_estimate")
  expect_true(is.na(res$estimate))
  expect_identical(res$interval_type, "one-sided")
  expect_equal(res$last_sighting, 100)
  expect_lt(abs(res$recovery_extension -
                  res$recovery_range * (0.05^(-1 / 5) - 1)), 1e-12)
})

test_that("marshall1997 validates its inputs", {
  sd <- dead_bones()
  expect_error(marshall1997(sd, alpha = 0), "must be in \\(0, 1\\)")
  expect_error(marshall1997(sd, recovery = 1), "must be a function")
  expect_error(marshall1997(sd, function(t) -t), "non-negative")
  expect_error(marshall1997(sd, function(t) 0 * t), "integrates to zero")
  expect_error(marshall1997(sd, max_time = 50), "max_time")
  expect_error(marshall1997(sd, test_year = 90), "later than the last sighting")
  one <- sighting_data(data.frame(time = 1, count = 1))
  expect_error(marshall1997(one), "at least 2")
})
