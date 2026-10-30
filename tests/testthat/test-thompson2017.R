# Alaotra Grebe inputs, Table 1 of Thompson et al. (2017)
grebe_records <- function() {
  data.frame(
    year = c(1929, 1960, 1963, 1969, 1970, 1971, 1972, 1982, 1985, 1986, 1988),
    pci_lower = c(0.99, 0.95, 0.75, 0.75, 0.10, 0.10, 0.60, 0.60, 0.20, 0.20, 0.20),
    pci_upper = c(1.00, 0.99, 0.94, 0.94, 0.40, 0.40, 0.80, 0.80, 0.80, 0.70, 0.50)
  )
}
grebe_surveys <- function() {
  y <- c(1989, 1990, 1993, 1994, 1997, 1998, 1999, 2000, 2004, 2009)
  data.frame(year = y,
             eps_lower = ifelse(y == 2000, 0.70, 0.80),
             eps_upper = ifelse(y == 2000, 0.90, 0.95),
             pi_lower = 0.90, pi_upper = 0.95, pr_lower = 0.70, pr_upper = 0.90)
}
grebe_passive <- c(eps_lower = 0, eps_upper = 0.05, pi_lower = 0.10,
                   pi_upper = 0.65, pr_lower = 0.40, pr_upper = 0.60)
# P(X0) in 0.95-0.99 (section 6)
grebe_fit <- function(records = grebe_records(), surveys = grebe_surveys(), ...) {
  thompson2017(records, surveys, grebe_passive, p_initial = c(0.95, 0.99),
               end_year = 2017, data_out = TRUE, ...)
}

test_that("thompson2017 iterates equation 5 with midpoints and interval bounds", {
  d <- grebe_fit(n_sim = 0)
  expect_equal(d$time, 1929:2017)
  at <- function(y, col = "p_extant") d[[col]][d$time == y]
  # 1929 record (eq. 1), then 30 passive years 1930-1959 (eqs. 3-4)
  p1929 <- 0.995 + (1 - 0.995) * 0.97
  expect_equal(at(1929), p1929)
  expect_equal(at(1959), p1929 * (1 - 0.025 * 0.375 * 0.5)^30)
  low1929 <- 0.99 + 0.01 * 0.95
  expect_equal(at(1959, "lower"), low1929 * (1 - 0.05 * 0.65 * 0.60)^30)
  expect_equal(at(1959, "upper"), 1)
  # first dedicated survey after the 1988 record
  expect_equal(at(1989), at(1988) * (1 - 0.875 * 0.925 * 0.8))
  expect_equal(at(1989, "upper"), at(1988, "upper") * (1 - 0.8 * 0.9 * 0.7))
  expect_true(all(d$lower <= d$p_extant & d$p_extant <= d$upper))
})

test_that("thompson2017 matches Fig. 1 of Thompson et al. (2017)", {
  d <- grebe_fit(n_sim = 0)
  at <- function(y, col = "p_extant") d[[col]][d$time == y]
  # values read from Fig. 1 (to about 0.01)
  expect_lt(abs(at(1959) - 0.87), 0.01)
  expect_lt(abs(at(1959, "lower") - 0.555), 0.01)
  expect_lt(abs(at(1981, "lower") - 0.83), 0.01)
  expect_lt(abs(at(1989) - 0.35), 0.01)
  expect_lt(abs(at(1990) - 0.12), 0.01)
  expect_lt(abs(at(1990, "upper") - 0.25), 0.01)
  expect_lt(abs(at(1994, "upper") - 0.06), 0.01)
  expect_lt(at(2000), 0.001)
})

test_that("thompson2017 matches Fig. 2 (hypothetical record in 2000)", {
  rec <- rbind(grebe_records(),
               data.frame(year = 2000, pci_lower = 0.2, pci_upper = 0.5))
  sv <- grebe_surveys()
  d <- grebe_fit(rec, sv[sv$year != 2000, ], n_sim = 0)
  at <- function(y, col = "p_extant") d[[col]][d$time == y]
  # values read from Fig. 2 (to about 0.01)
  expect_lt(abs(at(2001) - 0.35), 0.01)
  expect_lt(abs(at(2004) - 0.12), 0.01)
  expect_lt(abs(at(2004, "upper") - 0.25), 0.01)
  expect_lt(abs(at(2009) - 0.04), 0.01)
  expect_lt(abs(at(2017, "upper") - 0.125), 0.01)
})

test_that("thompson2017 Monte Carlo band matches Fig. 1", {
  set.seed(2017)
  d <- grebe_fit(n_sim = 20000)
  at <- function(y, col) d[[col]][d$time == y]
  # dark band read from Fig. 1 (to about 0.02)
  expect_lt(abs(at(1959, "mc_lower") - 0.835), 0.02)
  expect_lt(abs(at(1990, "mc_lower") - 0.07), 0.02)
  expect_lt(abs(at(1990, "mc_upper") - 0.18), 0.02)
  expect_true(all(d$mc_lower <= d$mc_upper))
  expect_true(all(d$mc_lower >= 0 & d$mc_upper <= 1))
})

test_that("thompson2017 passive-only projection follows equation 14", {
  # Section 8.2: eps' = 0.05, p'(r) = p'(i) = 0.8, P(X0) = 0.10
  d <- thompson2017(NULL, NULL, c(eps = 0.05, pi = 0.8, pr = 0.8),
                    p_initial = 0.1, start_year = 1, end_year = 30,
                    n_sim = 0, data_out = TRUE)
  alpha <- 0.05 * 0.8 * 0.8
  expect_equal(d$p_extant, (1 - alpha)^(1:30) * 0.1)
  # equation 15: t >= log(0.05 / 0.1) / log(1 - alpha) = 21.3 ("at least 21 years")
  expect_equal(min(d$time[d$p_extant <= 0.05]), 22)
  expect_equal(round(log(0.05 / 0.1) / log(1 - alpha), 1), 21.3)
})

test_that("thompson2017 summary object and edge cases", {
  res <- thompson2017(grebe_records(), grebe_surveys(), grebe_passive,
                      p_initial = c(0.95, 0.99), n_sim = 0)
  expect_s3_class(res, "ede_extant")
  expect_equal(res$end_year, 2009)
  expect_equal(res$p_extinct, 1 - res$p_extant)
  expect_true(is.na(res$mc_lower))
  out <- capture.output(print.ede_extant(res))
  expect_match(out[2], "end of 2009")
  # no passive surveillance and p(ci) = 1 records
  d <- thompson2017(data.frame(year = 2000, pci = 1), NULL,
                    c(eps = 0, pi = 0, pr = 0), p_initial = 0.3,
                    end_year = 2005, n_sim = 0, data_out = TRUE)
  expect_equal(d$p_extant, rep(1, 6))
  expect_equal(d$type, c("record", rep("passive", 5)))
  # a complete, perfect survey drives the probability to zero (eq. 4)
  d <- thompson2017(NULL, data.frame(year = 1, eps = 1, pi = 1, pr = 1),
                    c(eps = 0, pi = 0, pr = 0), p_initial = 0.9,
                    n_sim = 0, data_out = TRUE)
  expect_equal(d$p_extant, 0)
})

test_that("thompson2017 validates its inputs", {
  rec <- grebe_records()
  sv <- grebe_surveys()
  expect_error(thompson2017(rec, sv, p_initial = 0.9), "passive")
  expect_error(thompson2017(rec, sv, grebe_passive), "p_initial")
  expect_error(thompson2017(rec, sv, grebe_passive, c(0.99, 0.95)), "lower, upper")
  expect_error(thompson2017(rec, sv, grebe_passive, 1.2), "\\[0, 1\\]")
  bad <- rec
  bad$pci_lower[1] <- 1.5
  expect_error(thompson2017(bad, sv, grebe_passive, 0.9), "pci")
  bad <- rec
  bad$pci_upper[2] <- 0.5
  expect_error(thompson2017(bad, sv, grebe_passive, 0.9), "exceed")
  expect_error(thompson2017(rec, sv[, -2], grebe_passive, 0.9), "eps")
  both <- rbind(rec, data.frame(year = 1989, pci_lower = 0.1, pci_upper = 0.2))
  expect_error(thompson2017(both, sv, grebe_passive, 0.9), "both a record")
  expect_error(thompson2017(rbind(rec, rec[1, ]), sv, grebe_passive, 0.9),
               "duplicates")
  expect_error(thompson2017(NULL, NULL, grebe_passive, 0.9), "start_year")
  expect_error(thompson2017(rec, sv, grebe_passive, 0.9, end_year = 1900),
               "earlier")
  expect_error(thompson2017(rec, sv, grebe_passive, 0.9, n_sim = -1), "n_sim")
})
