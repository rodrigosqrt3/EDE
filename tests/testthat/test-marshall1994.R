# Metrarabdotos n. sp. 5 (Marshall 1994, p. 462-463): 11 horizons with ages
# binned to 0.05 Myr, giving 10 gaps in stratigraphic order
# 0.00, 0.05, 0.00, 0.00, 1.10, 0.10, 0.00, 0.10, 0.00, 0.05 Myr.
# The paper prints the 8th gap as 0.05; it must be 0.10 to match both the
# ranked list (two 0.05 and two 0.10 gaps) and T = 14 (see test below).
# Positions are coded as Myr along the range; a row with count k is k
# horizons in the same 0.05-Myr bin.
metrarabdotos5 <- function() {
  sighting_data(data.frame(time = c(0, 0.05, 1.15, 1.25, 1.35, 1.40),
                           count = c(2, 3, 1, 2, 2, 1)))
}

test_that("gap ranks reproduce Table 1 (Marshall 1994)", {
  # rows N: l-u for C = 0.5, 0.8, 0.95 at confidence probabilities 0.95 and
  # 0.99, NA = no bound. Taken from Table 1; entries where the printed value
  # differs from equations 1-2 (N = 24, 39, 65, 70 for C = 0.5; N = 90, 95
  # for C = 0.95 at 0.99, OCR-damaged) are omitted.
  table1 <- list(
    `6`  = list(c(1, 6),   NULL,      c(3, NA),  c(2, NA),  c(4, NA),  c(4, NA)),
    `10` = list(c(2, 9),   c(1, 10),  c(5, NA),  c(4, NA),  c(8, NA),  c(7, NA)),
    `17` = list(c(5, 13),  c(3, 15),  c(10, 17), c(9, NA),  c(14, NA), c(13, NA)),
    `20` = list(c(6, 15),  c(4, 17),  c(12, 20), c(11, NA), c(17, NA), c(16, NA)),
    `25` = list(c(8, 18),  c(6, 20),  c(16, 25), c(14, 25), c(21, NA), c(20, NA)),
    `30` = list(c(10, 21), c(8, 23),  c(19, 29), c(18, 30), c(26, NA), c(25, NA)),
    `41` = list(c(14, 28), c(12, 30), c(28, 38), c(26, 40), c(36, NA), c(35, NA)),
    `50` = list(c(18, 33), c(16, 35), c(34, 46), c(32, 48), c(44, NA), c(43, NA)),
    `60` = list(c(22, 39), c(20, 41), c(42, 55), c(40, 56), c(53, NA), c(52, NA)),
    `75` = list(c(29, 47), c(26, 50), c(53, 67), c(51, 69), c(67, 75), c(66, NA)),
    `80` = list(c(31, 50), c(29, 52), c(57, 72), c(54, 74), c(72, 80), c(70, NA)),
    `100` = list(c(40, 61), c(37, 64), c(72, 89), c(69, 91), c(90, 100), c(89, NA))
  )
  settings <- expand.grid(cp = c(0.95, 0.99), C = c(0.5, 0.8, 0.95))
  for (n in names(table1)) {
    for (k in seq_len(nrow(settings))) {
      want <- table1[[n]][[k]]
      if (is.null(want)) next
      got <- marshall1994_ranks(as.integer(n), settings$C[k],
                                (1 - settings$cp[k]) / 2)
      expect_equal(unname(got), as.integer(want),
                   info = paste("N =", n, "C =", settings$C[k], "cp =", settings$cp[k]))
    }
  }
  # N = 6, C = 0.5 at 0.99 has no bounds; N = 5 at 0.95 neither (p. 463:
  # at least seven horizons for a 50% interval)
  expect_true(all(is.na(marshall1994_ranks(6L, 0.5, 0.005))))
  expect_true(all(is.na(marshall1994_ranks(5L, 0.5, 0.025))))
})

test_that("an upper limit on the 95% interval needs 72 gaps (p. 462, 464)", {
  has_upper <- vapply(1:120, function(n) {
    !is.na(marshall1994_ranks(n, 0.95, 0.025)[["upper"]])
  }, logical(1))
  expect_equal(which(has_upper)[1], 72L)
})

test_that("the binomial sums reproduce Table 2", {
  # Table 2: cumulative Binomial(10, C) probabilities, x = 0..10
  c50 <- c(0.001, 0.011, 0.055, 0.172, 0.377, 0.623, 0.828, 0.945, 0.989, 0.999, 1)
  c95 <- c(0, 0, 0, 0, 0, 0, 0.001, 0.012, 0.086, 0.401, 1)
  expect_equal(round(stats::pbinom(0:10, 10, 0.5), 3), c50)
  expect_equal(round(stats::pbinom(0:10, 10, 0.95), 3), c95)
})

test_that("marshall1994 reproduces the Metrarabdotos n. sp. 5 example", {
  sd <- metrarabdotos5()
  # 50% interval: between the 2nd and 9th smallest gaps, 0.00-0.10 Myr
  res <- marshall1994(sd, alpha = 0.5, horizons = "count")
  expect_equal(unname(res$gap_ranks), c(2L, 9L))
  expect_equal(unname(res$ci_length), c(0, 0.10))
  expect_equal(res$upper, 1.40 + 0.10)
  expect_equal(res$upper_min, 1.40)
  # 95% interval: lower limit the 8th smallest gap (0.10 Myr), no upper limit
  res <- marshall1994(sd, alpha = 0.05, horizons = "count")
  expect_equal(unname(res$gap_ranks), c(8L, NA))
  expect_equal(res$ci_length[["lower"]], 0.10)
  expect_equal(res$upper, Inf)
  # 80% interval (Table 5): 0.00 to more than the largest gap
  res <- marshall1994(sd, alpha = 0.2, horizons = "count")
  expect_equal(res$ci_length[["lower"]], 0)
  expect_equal(res$upper, Inf)
  expect_equal(res$n_gaps, 10L)
})

test_that("marshall1994_trend reproduces the transposition example (p. 463)", {
  tr <- marshall1994_trend(metrarabdotos5(), horizons = "count")
  expect_equal(tr$gaps, c(0, 0.05, 0, 0, 1.10, 0.10, 0, 0.10, 0, 0.05),
               tolerance = 1e-12)
  # 14 transpositions + 5 * 4 / 4 + 1/2 + 1/2 = 20, T_max / 2 = 22.5
  expect_equal(tr$transpositions, 20)
  expect_equal(tr$t_max / 2, 22.5)
  expect_identical(tr$direction, "increasing")
  expect_gt(tr$p_value, 0.2)
})

test_that("the transposition p-value reproduces Table 3 critical values", {
  # Table 3: critical T at p = 0.01, 0.05, 0.10, 0.20 for N = 3..20
  table3 <- rbind(
    `3` = c(NA, NA, NA, NA), `4` = c(NA, NA, 0, 0), `5` = c(NA, 0, 1, 1),
    `6` = c(0, 1, 2, 3), `7` = c(1, 3, 4, 5), `8` = c(3, 5, 6, 8),
    `9` = c(5, 8, 9, 11), `10` = c(8, 11, 12, 14), `11` = c(11, 14, 16, 18),
    `12` = c(14, 18, 20, 23), `13` = c(18, 22, 24, 27), `14` = c(22, 27, 29, 33),
    `15` = c(26, 32, 35, 38), `16` = c(31, 37, 41, 45), `17` = c(36, 43, 47, 51),
    `18` = c(42, 50, 54, 58), `19` = c(49, 56, 61, 66), `20` = c(54, 64, 69, 74)
  )
  # 67 of the 72 entries match exactly. Five differ by one at borderline
  # probabilities (exact two-sided p at the tabulated value / at ours):
  # N = 13, p = .01: 0.0101 / 0.0067; N = 13, p = .10: 0.0763 / 0.09998;
  # N = 19, p = .01: 0.0103 / 0.0083; N = 19, p = .05: 0.0412 / 0.0491;
  # N = 20, p = .01: 0.0074 / 0.0091 (the table is conservative there).
  off_by_one <- c("13 0.01", "13 0.1", "19 0.01", "19 0.05", "20 0.01")
  levels <- c(0.01, 0.05, 0.10, 0.20)
  for (n in rownames(table3)) {
    p_at <- pmin(1, 2 * cumsum(marshall1994_inversions(as.integer(n))))
    for (k in seq_along(levels)) {
      crit <- if (any(p_at <= levels[k])) max(which(p_at <= levels[k])) - 1 else NA_real_
      key <- paste(n, levels[k])
      if (key %in% off_by_one) {
        expect_equal(abs(crit - unname(table3[n, k])), 1, info = key)
      } else {
        expect_equal(crit, unname(table3[n, k]), info = key)
      }
    }
  }
})

test_that("marshall1994_trend detects monotone gap sequences", {
  inc <- sighting_data(data.frame(time = cumsum(1:12), count = 1))
  tr <- marshall1994_trend(inc)
  expect_equal(tr$transpositions, 0)
  expect_identical(tr$direction, "increasing")
  expect_lt(tr$p_value, 0.01)
  dec <- sighting_data(data.frame(time = -cumsum(1:12), count = 1))
  tr <- marshall1994_trend(dec)
  expect_equal(tr$transpositions, tr$t_max)
  expect_identical(tr$direction, "decreasing")
  expect_lt(tr$p_value, 0.01)
})

test_that("marshall1994 uses distinct horizons by default and validates input", {
  sd <- metrarabdotos5()
  expect_equal(marshall1994(sd)$n_gaps, 5L)
  res <- marshall1994(sd, alpha = 0.5)
  expect_s3_class(res, "ede_estimate")
  expect_true(is.na(res$estimate))
  expect_identical(res$interval_type, "one-sided")
  expect_error(marshall1994(sd, alpha = 0), "must be in \\(0, 1\\)")
  expect_error(marshall1994(sd, conf_prob = 1), "conf_prob")
  expect_error(marshall1994(sd, horizons = "x"))
  one <- sighting_data(data.frame(time = 1, count = 1))
  expect_error(marshall1994(one), "at least 2 find horizons")
  expect_error(marshall1994_trend(sighting_data(data.frame(time = 1:2, count = 1))),
               "at least 2 gaps")
})
