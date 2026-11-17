rec <- function(years) sighting_data(data.frame(year = years, sightings = 1))

eskimo <- c(1946, 1947, 1948, 1950, 1955, 1956, 1959, 1960, 1961, 1962, 1963,
            1964, 1968, 1970, 1972, 1973, 1974, 1976, 1977, 1980, 1981, 1982,
            1983, 1985)

# Table 2 of Vogel et al. (2009): n, r and the published p-value of the PPCC
# test of the uniform distribution for 27 bird records
vg_n <- c(6, 8, 17, 8, 24, 15, 28, 9, 13, 19, 7, 11, 19, 10, 18, 11, 7, 6, 11,
          11, 9, 9, 6, 10, 8, 12, 26)
vg_r <- c(.984, .994, .986, .917, .996, .979, 1, .991, .992, .993, .988, .98,
          .997, .996, .969, .988, .947, .95, .995, .9, .749, .875, .911, .933,
          .979, .978, .992)
vg_p <- c(.921, .999, .651, .08, .905, .432, 1, .967, .934, .893, .951, .683,
          .996, .999, .154, .863, .287, .383, .999, .017, .001, .015, .128,
          .086, .787, .581, .745)

test_that("vogel2009 reproduces the PPCC of the Eskimo curlew record", {
  res <- vogel2009(rec(eskimo), nsim = 2e4)
  expect_s3_class(res, "htest")
  expect_equal(unname(res$statistic), 0.996, tolerance = 0.001)
  expect_equal(unname(res$parameter), 24)
  expect_gt(res$p.value, 0.5)   # no evidence against the uniform model
})

test_that("null distribution agrees with Vogel & Kroll (1989) and Table 2", {
  # 95% of samples of length 24 have r above about 0.96 (Vogel et al., p. 711)
  expect_equal(unname(stats::quantile(EDE:::ppcc_uniform_null(24L, 5e4), 0.05)),
               0.96, tolerance = 0.01)
  # published p-values agree in the lower tail, where rejections are decided;
  # in the upper tail they differ because r is given to three decimals
  mine <- mapply(function(r, n) EDE:::ppcc_uniform_p(r, n, 5e4), vg_r, vg_n)
  low <- vg_p < 0.2
  expect_lt(max(abs(mine[low] - vg_p[low])), 0.02)
  # the three records rejected at 5% are the same
  expect_equal(mine < 0.05, vg_p < 0.05)
})

test_that("field significance reproduces Vogel et al. (2009, Section 6)", {
  res <- ede_field_test(vg_p)
  expect_equal(res$rejections[res$test == "Bonferroni-type"], 1)
  expect_equal(res$rejections[grepl("discovery", res$test)], 1)
  bin <- res[res$test == "Binomial", ]
  expect_equal(bin$rejections, 3)
  expect_equal(bin$p_value, 0.15, tolerance = 0.01)
  expect_equal(ede_field_test(vg_p, alpha = 0.01)$p_value[3], 0.24,
               tolerance = 0.01)
  # the paper reads this value from Figure 4, hence the looser tolerance
  expect_equal(ede_field_test(vg_p, alpha = 0.1)$p_value[3], 0.13,
               tolerance = 0.05)
  ppcc <- res[grepl("PPCC", res$test), ]
  expect_equal(ppcc$statistic, 0.9597, tolerance = 1e-3)
  expect_equal(ppcc$p_value, 0.022, tolerance = 0.15)
})

test_that("statistic is invariant to the origin and unit of time", {
  a <- vogel2009(rec(eskimo), nsim = 1e4)
  b <- vogel2009(rec(5 * eskimo - 2000), nsim = 1e4)
  expect_equal(a$statistic, b$statistic)
  expect_equal(a$p.value, b$p.value)
})

test_that("size is controlled and a non-uniform record is rejected", {
  skip_on_cran()
  set.seed(41)
  p0 <- replicate(300, vogel2009(rec(runif(20)), nsim = 1e4)$p.value)
  expect_lt(mean(p0 < 0.05), 0.10)
  p1 <- replicate(100, vogel2009(rec(rexp(20)), nsim = 1e4)$p.value)
  expect_gt(mean(p1 < 0.05), 0.2)
})

test_that("does not touch the global random number stream", {
  set.seed(1)
  before <- .Random.seed
  vogel2009(rec(eskimo), nsim = 1e4)
  expect_identical(.Random.seed, before)
})

test_that("input validation", {
  expect_error(vogel2009(1), "sighting_data")
  expect_error(vogel2009(rec(c(1, 2, 3))), "at least 4")
  expect_error(vogel2009(rec(eskimo), nsim = 10), "nsim")
  expect_error(ede_field_test(c(0.1, 0.2)), "at least 4")
  expect_error(ede_field_test(c(0.1, 0.2, 0.3, 2)), "\\[0, 1\\]")
  expect_error(ede_field_test(vg_p, alpha = 0), "alpha")
})

test_that("sighting times that are all equal are refused", {
  flat <- sighting_data(data.frame(time = 3, count = 8))
  expect_error(vogel2009(flat, events = TRUE, nsim = 2000),
               "must not all be equal")
})

test_that("the null restores the global random stream", {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  EDE:::ppcc_uniform_null(7, 555)
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})
