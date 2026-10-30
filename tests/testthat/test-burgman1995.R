seal_record <- function() {
  sighting_data(data.frame(
    year = c(1915, 1922, 1932, 1948, 1952),
    sightings = 1
  ))
}

test_that("burgman1995 reproduces the monk seal values in Table 2", {
  # The first sighting (1915) defines the origin: CT = 77, Ce = 37, N = 4.
  eq2 <- burgman1995(seal_record(), test_year = 1992, data_out = TRUE)
  expect_equal(eq2$chance[eq2$time == 1992], (37 / 77)^4)
  expect_equal(round(eq2$chance[eq2$time == 1992], 3), 0.053)

  eq4 <- burgman1995(seal_record(), test_year = 1992, data_out = TRUE,
                     equation = 4)
  expect_equal(round(eq4$chance[eq4$time == 1992], 3), 0.258)
})

test_that("burgman1995 equation 4 reproduces the Fig. 1b example", {
  # CT = 16 periods, N = 10 observations, longest empty run r = 9
  expect_equal(round(EDE:::burgman_equation4_pvalue(ct = 16, n = 10, r = 9), 5),
               0.00167)
})

test_that("burgman1995 returns the first rejected period as an upper bound", {
  res <- burgman1995(seal_record(), alpha = 0.05, test_year = 2100)
  expect_true(is.na(res$estimate))
  expect_equal(res$interval_type, "one-sided")
  # (37 / (T - 1915))^4 <= 0.05 first holds at T = 1994
  expect_equal(res$upper, 1994)
})

test_that("burgman1995 equation 2 uses only the final run of empty periods", {
  d <- sighting_data(data.frame(year = c(1900:1902, 1950:1952), sightings = 1))
  eq2 <- burgman1995(d, test_year = 2100)
  eq4 <- burgman1995(d, test_year = 2100, equation = 4)
  # the long gap in the middle of the record drives equation 4 only
  expect_equal(eq4$upper, 1953)
  expect_gt(eq2$upper, 1990)
})

test_that("burgman1995 validates its inputs", {
  d <- data.frame(years = c(1900, 1910), sightings = c(1, 1))
  sd <- sighting_data(d)
  expect_error(burgman1995(sd, alpha = 0, test_year = 2000), "must be in \\(0, 1\\)")
  expect_error(burgman1995(sd, alpha = 1, test_year = 2000), "must be in \\(0, 1\\)")
  expect_error(burgman1995(sd, alpha = 0.05), "test_year")
  expect_error(burgman1995(sd, alpha = 0.05, test_year = "2000"), "test_year")
  expect_error(burgman1995(sd, alpha = 0.05, test_year = c(2000, 2001)), "test_year")
  expect_error(burgman1995(sd, test_year = 2000, equation = 5), "2, 3 or 4")
  expect_error(burgman1995(sd, test_year = 2000, equation = c(2, 3)), "2, 3 or 4")

  fractional <- sighting_data(data.frame(time = c(1900.5, 1910), count = 1))
  expect_error(burgman1995(fractional, test_year = 2000), "whole-number")
})

test_that("burgman1995 rejects test_year equal to the last sighting", {
  d <- data.frame(
    years = c(1907, 1910, 1915, 1916, 1920, 1925, 1930, 1931),
    sightings = c(1, 1, 3, 4, 3, 1, 2, 1)
  )
  expect_error(burgman1995(sighting_data(d), 0.05, test_year = 1931), "later than")
})

test_that("burgman1995 data_out = TRUE returns the full p-value curve", {
  d <- data.frame(
    years = c(1907, 1910, 1915, 1916, 1920, 1925, 1930, 1931),
    sightings = c(1, 1, 3, 4, 3, 1, 2, 1)
  )
  curve <- burgman1995(sighting_data(d), alpha = 0.05, test_year = 1935, data_out = TRUE)
  expect_named(curve, c("time", "chance"))
  expect_equal(nrow(curve), 4L)
})

test_that("burgman1995 warns and returns NA when alpha is not reached", {
  expect_warning(
    res <- burgman1995(seal_record(), alpha = 0.05, test_year = 1960),
    "never falls"
  )
  expect_true(is.na(res$upper))
})

test_that("burgman equation 4 stays finite for large counts", {
  p <- EDE:::burgman_equation4_pvalue(ct = 100, n = 200, r = 10)
  expect_true(is.finite(p))
  expect_gte(p, 0)
  expect_lt(p, 1e-4)
  # agrees with a direct simulation of 60 records in 100 periods (~0.078)
  expect_lt(abs(EDE:::burgman_equation4_pvalue(ct = 100, n = 60, r = 10) - 0.078), 3e-3)
})

test_that("burgman internal helper functions handle edge cases", {
  expect_equal(EDE:::burgman_equation4_pvalue(ct = 10, n = 0, r = 2), 1.0)
  expect_equal(EDE:::burgman_equation4_pvalue(ct = 10, n = 5, r = 0), 1.0)
  expect_equal(EDE:::burgman_equation4_pvalue(ct = 10, n = 5, r = 10), 1.0)
  expect_equal(EDE:::burgman_equation4_pvalue(ct = 10, n = 5, r = 7), 0.00876, tolerance = 1e-5)
  d_zero <- data.frame(time = 1900:1905, count = rep(0, 6))
  expect_equal(EDE:::burgman_chance(d_zero), 1.0)
  expect_equal(EDE:::burgman_equation2_chance(d_zero), 1.0)
  expect_equal(EDE:::burgman_equation3_chance(d_zero), 1.0)
})

test_that("burgman1995 equation 3 reproduces the Fig. 1b example", {
  # Burgman et al. (1995), Fig. 1b: CT = 16, n1 = 5 occupied periods, r = 9
  expect_equal(round(EDE:::burgman_equation3_pvalue(ct = 16, n1 = 5, r = 9), 5),
               0.02885)
})

test_that("burgman1995 equation 3 reproduces Table 2 for the ferret", {
  # Table 1: yearly ferret counts 1972-1984, last period 1992. Observation
  # began in January 1972, so a zero-count 1971 row makes 1972 a period:
  # CT = 21, n1 = 11, r = 8. Table 2: eq. 2 0.000, eq. 3 0.003, eq. 4 0.000.
  ferret <- sighting_data(data.frame(
    year = c(1971:1977, 1979, 1981:1984),
    count = c(0, 4, 6, 2, 3, 3, 1, 1, 2, 3, 1, 2)
  ))
  at_1992 <- function(eq) {
    curve <- burgman1995(ferret, test_year = 1992, data_out = TRUE, equation = eq)
    round(curve$chance[curve$time == 1992], 3)
  }
  expect_equal(at_1992(3), 0.003)
  expect_equal(at_1992(2), 0)
  expect_equal(at_1992(4), 0)
})

test_that("burgman1995 equation 3 reproduces Table 2 for the stick-nest rats", {
  # Table 2: L. conditor eq. 3 0.139 (eq. 4 0.067), L. apicalis eq. 3 0.350
  # (eq. 4 0.268). The paper's CT and r are not recoverable exactly from the
  # dates in Table 1 (which give CT = 153); the values below are those for
  # which the verified equation 4 also reproduces the paper.
  expect_equal(round(EDE:::burgman_equation3_pvalue(155, 12, 47), 3), 0.139)
  expect_equal(round(EDE:::burgman_equation4_pvalue(155, 15, 47), 3), 0.067)
  expect_equal(round(EDE:::burgman_equation3_pvalue(156, 7, 55), 3), 0.350)
  expect_equal(round(EDE:::burgman_equation4_pvalue(156, 8, 55), 3), 0.268)
})

test_that("burgman1995 equation 3 for the monk seal (paper value not reproduced)", {
  # CT = 77, n1 = 4, r = 40 reproduce the paper's eq. 2 (0.053) and eq. 4
  # (0.258) values but give 0.244 for eq. 3; Table 2 lists 0.258.
  eq3 <- burgman1995(seal_record(), test_year = 1992, data_out = TRUE, equation = 3)
  expect_equal(round(eq3$chance[eq3$time == 1992], 3), 0.244)
})

test_that("burgman1995 equation 3 ignores counts and agrees with enumeration", {
  d1 <- sighting_data(data.frame(year = c(1899, 1900, 1905, 1906), count = c(0, 1, 3, 1)))
  d2 <- sighting_data(data.frame(year = c(1899, 1900, 1905, 1906), count = c(0, 1, 1, 1)))
  c1 <- burgman1995(d1, test_year = 1915, data_out = TRUE, equation = 3)
  c2 <- burgman1995(d2, test_year = 1915, data_out = TRUE, equation = 3)
  expect_equal(c1$chance, c2$chance)

  # brute force over all placements of n1 occupied cells among ct
  brute <- function(ct, n1, r) {
    occ <- utils::combn(ct, n1)
    mean(apply(occ, 2, function(o) {
      x <- rep(0, ct); x[o] <- 1
      EDE:::burgman_longest_empty_run(x) >= r
    }))
  }
  for (x in list(c(10, 3, 3), c(12, 4, 2), c(9, 2, 5), c(8, 1, 7))) {
    expect_lt(abs(EDE:::burgman_equation3_pvalue(x[1], x[2], x[3]) -
                    brute(x[1], x[2], x[3])), 1e-12)
  }
  expect_equal(EDE:::burgman_equation3_pvalue(10, 3, 0), 1)
  expect_equal(EDE:::burgman_equation3_pvalue(10, 3, 8), 0)
  expect_true(is.finite(EDE:::burgman_equation3_pvalue(2000, 300, 40)))
})

test_that("burgman1995 reproduces Burgman et al. (2000) Table 2.1", {
  # Western Australian Acacia: the first collection year is the first period
  # (zero-count row for the year before) and 1993 is the last period.
  # Columns 'Solow' (eq. 2) and 'Grimson' (eq. 4) of Table 2.1.
  acacia <- list(
    ancistrophylla = list(y = c(37, 39, 68, 70, 71, 71, 73, 76, 79, 89), p = c(.483, .011)),
    excentrica = list(y = c(78, 83), p = c(.141, .398)),
    gemina = list(y = c(38, 64, 70, 71, 71, 80, 80, 82, 82, 82, 87, 91), p = c(.646, .009)),
    kingiana = list(y = 23, p = c(.014, .028)),
    megacephala = list(y = c(32, 56, 61, 63, 64, 64, 67, 67, 70, 70, 70, 70, 71, 71, 72, 76),
                       p = c(.006, .009)),
    prismifolia = list(y = c(1, 33), p = c(.126, .374)),
    subsessilis = list(y = c(31, 66, 71, 76, 82, 87), p = c(.549, .062))
  )
  for (sp in acacia) {
    tb <- table(1900 + sp$y)
    yrs <- as.numeric(names(tb))
    d <- sighting_data(data.frame(time = c(min(yrs) - 1, yrs),
                                  count = c(0, as.numeric(tb))))
    for (i in 1:2) {
      curve <- burgman1995(d, test_year = 1993, data_out = TRUE, equation = c(2, 4)[i])
      expect_equal(round(curve$chance[curve$time == 1993], 3), sp$p[i])
    }
  }
})
