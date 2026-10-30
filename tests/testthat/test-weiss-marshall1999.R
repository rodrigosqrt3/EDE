# Weiss & Marshall (1999), Table 1: presences of 16 benthic foraminifera at 36
# sampling points of the Cenomanian-Turonian at Eastbourne (heights in m above
# the first sample). Each vector gives the sampling-point indices j with a
# presence; every other sampling point is an absence.
wm_heights <- c(0, 1.4, 1.6, 2.3, 3.9, 4.1, 4.7, 5.1, 5.4, 5.9, 6.1, 6.6, 7,
                7.4, 7.9, 8.6, 9.2, 9.9, 10.2, 10.4, 10.8, 11.2, 11.4, 11.7,
                11.9, 12.5, 14.1, 14.9, 15.4, 16.3, 18.4, 20.7, 22.5, 25.8,
                27.6, 31.4)
wm_presence <- list(
  `20` = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
  `21` = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
  `22` = c(1, 2, 3, 4, 5, 6, 7, 9, 10),
  `23` = c(1, 2, 3, 4, 5, 6, 7, 8, 10),
  `24` = c(1, 2, 3, 6, 7, 8, 9, 10),
  `25` = c(11, 12, 14),
  `26` = c(2, 13, 16),
  `27` = c(4, 9, 10, 11, 17),
  `28` = c(15, 17),
  `29` = c(2, 4, 15, 18),
  `30` = c(6, 18),
  `43` = c(1, 2, 4, 5, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 23, 24),
  `44` = c(1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
           22, 23, 24, 25),
  `45` = c(1, 2, 3, 4, 5, 7, 10, 11, 12, 13, 15, 16, 17, 18, 19, 20, 22, 23, 25),
  `46` = c(1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
           21, 22, 23, 24, 25, 27),
  `47` = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 15, 16, 19, 21, 22, 23, 25,
           26, 27, 29)
)
wm_record <- function(sp) {
  count <- rep(0, length(wm_heights))
  count[wm_presence[[sp]]] <- 1
  sighting_data(data.frame(time = wm_heights, count = count))
}
wm_groups <- list(as.character(20:24), as.character(25:30), as.character(43:47))
wm_prior <- function(sp) {
  g <- wm_groups[[which(vapply(wm_groups, function(x) sp %in% x, logical(1)))]]
  weiss_marshall1999_prior(lapply(setdiff(g, sp), wm_record))
}
# residual summaries as tabulated in the paper: mean, SD, S.5, S.95 of S - t[m]
wm_summary <- function(res) {
  c(res$posterior_mean - res$lower, res$posterior_sd,
    res$estimate - res$lower, res$upper - res$lower, res$p_mean, res$p_sd)
}
expect_matches_table <- function(got, want) {
  # the paper prints 3 significant figures
  expect_true(all(abs(got - want) <= 0.006 * abs(want) + 1e-3),
              info = paste(signif(got, 4), collapse = " "))
}

test_that("weiss_marshall1999_prior reproduces the pooled priors of Table 2", {
  # Table 2, columns c and d
  table2_cd <- rbind(
    `20` = c(31, 4), `21` = c(31, 4), `22` = c(31, 3), `23` = c(31, 3),
    `24` = c(32, 2), `25` = c(11, 70), `26` = c(11, 68), `27` = c(9, 69),
    `28` = c(12, 66), `29` = c(10, 67), `30` = c(12, 65), `43` = c(85, 17),
    `44` = c(81, 20), `45` = c(85, 16), `46` = c(79, 20), `47` = c(82, 15)
  )
  for (sp in rownames(table2_cd)) {
    expect_equal(unname(wm_prior(sp)), table2_cd[sp, ], info = sp)
  }
})

test_that("weiss_marshall1999 reproduces Table 2 (l = 800 m)", {
  # Weiss & Marshall (1999), Table 2: posterior mean, SD, S.5, S.95 of the
  # distance beyond the last presence, posterior mean and SD of p
  table2 <- rbind(
    `20` = c(.270, .175, .262, .498, .908, .043),
    `22` = c(.180, .194, .125, .614, .905, .045),
    `24` = c(.180, .194, .125, .614, .905, .045),
    `25` = c(730, 797, 481, 2320, .122, .031),
    `26` = c(748, 798, 501, 2340, .123, .031),
    `27` = c(757, 799, 510, 2350, .124, .031),
    `29` = c(765, 799, 519, 2360, .125, .031),
    `43` = c(.482, 4.30, .180, 2.01, .821, .034),
    `44` = c(.738, 5.51, .454, 2.11, .823, .034),
    `46` = c(.896, 23.2, .466, 1.29, .824, .034),
    `47` = c(5.85, 86.6, .681, 3.86, .822, .034)
  )
  for (sp in rownames(table2)) {
    res <- weiss_marshall1999(wm_record(sp), wm_prior(sp), mean_longevity = 800)
    expect_matches_table(wm_summary(res), table2[sp, ])
  }
})

test_that("weiss_marshall1999 reproduces the sensitivity analysis of Table 3", {
  # Table 3, first block: l = 40 and l = 2
  res <- weiss_marshall1999(wm_record("25"), c(11, 70), mean_longevity = 40)
  expect_matches_table(wm_summary(res), c(16.5, 30.0, 3.70, 79.0, .134, .035))
  res <- weiss_marshall1999(wm_record("25"), c(11, 70), mean_longevity = 2)
  expect_matches_table(wm_summary(res), c(1.25, 1.25, .899, 3.56, .145, .036))
  res <- weiss_marshall1999(wm_record("43"), c(85, 17), mean_longevity = 2)
  expect_matches_table(wm_summary(res), c(.327, .451, .154, 1.33, .823, .034))
  # changing the prior for p of species 20
  res <- weiss_marshall1999(wm_record("20"), c(11, 24), mean_longevity = 800)
  expect_matches_table(wm_summary(res), c(4.77, 79.7, .476, 2.75, .441, .076))
  # flat Beta(1, 1) prior, last block
  res <- weiss_marshall1999(wm_record("20"), c(1, 1), mean_longevity = 800)
  expect_matches_table(wm_summary(res), c(.281, 2.37, .264, .516, .904, .088))
  res <- weiss_marshall1999(wm_record("43"), c(1, 1), mean_longevity = 800)
  expect_matches_table(wm_summary(res), c(5.36, 86.4, .323, 2.40, .743, .088))
})

test_that("weiss_marshall1999 returns a one-sided ede_estimate on the time scale", {
  res <- weiss_marshall1999(wm_record("20"), c(31, 4), mean_longevity = 800,
                            alpha = 0.05)
  expect_s3_class(res, "ede_estimate")
  expect_identical(res$interval_type, "one-sided")
  expect_equal(res$lower, 5.4)
  expect_true(res$estimate > 5.4 && res$estimate < res$upper)
  # the 95% bound sits just below the next sampling point (paper, p. 445)
  expect_true(res$upper < 5.9)
})

test_that("persistence_prob is the posterior upper tail", {
  sd <- wm_record("43")
  res <- weiss_marshall1999(sd, c(85, 17), mean_longevity = 800, test_year = 11.7 + 2.007)
  expect_lt(abs(res$persistence_prob - 0.05), 2e-3)
  res <- weiss_marshall1999(sd, c(85, 17), 800, test_year = res$estimate)
  expect_lt(abs(res$persistence_prob - 0.5), 1e-8)
  expect_error(weiss_marshall1999(sd, c(85, 17), 800, test_year = 11.7),
               "later than the last sighting")
})

test_that("weiss_marshall1999 handles a presence at the last sample", {
  sd <- sighting_data(data.frame(time = 0:4, count = c(1, 0, 1, 1, 1)))
  res <- weiss_marshall1999(sd, c(2, 2), mean_longevity = 10)
  # beyond the last sample the posterior is the exponential prior
  expect_equal(res$estimate, 4 + 10 * log(2))
  expect_equal(res$posterior_mean, 4 + 10)
})

test_that("weiss_marshall1999 validates its inputs", {
  sd <- wm_record("20")
  expect_error(weiss_marshall1999(sd, c(31, 0), 800), "prior_p")
  expect_error(weiss_marshall1999(sd, 31, 800), "prior_p")
  expect_error(weiss_marshall1999(sd, c(31, 4), -1), "mean_longevity")
  expect_error(weiss_marshall1999(sd, c(31, 4), Inf), "mean_longevity")
  expect_error(weiss_marshall1999(sd, c(31, 4), 800, alpha = 1), "must be in \\(0, 1\\)")
  empty <- sighting_data(data.frame(time = 0:3, count = 0))
  expect_error(weiss_marshall1999(empty, c(1, 1), 10), "no presences")
  expect_error(weiss_marshall1999_prior(list()), "non-empty")
  expect_error(weiss_marshall1999_prior(list(empty)), "at least one presence")
})
