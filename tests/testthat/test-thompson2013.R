# Caribbean monk seal (Thompson et al. 2013, Section 3.3): S = 5, T_N = 38,
# T = 97 (1915-2011).
seal13 <- function() {
  sighting_data(data.frame(year = c(1915, 1922, 1932, 1948, 1952), sightings = 1))
}
th_seal <- function(...) thompson2013(seal13(), test_year = 2011, ...)

test_that("thompson2013 reproduces the deterministic monk seal estimates", {
  # eq 18 with p = 1 - 5/38: 0.0012; p = 0.92: 0.0219 (p. 19)
  res <- th_seal(p = 1 - 5 / 38)
  expect_equal(round(res$posterior, 4), 0.0012)
  expect_equal(round(th_seal(p = 0.92)$posterior, 4), 0.0219)
  expect_equal(res$average, "deterministic")
  expect_equal(res$prior, 38 / 97)
  # eq 18 in closed form
  expect_equal(res$posterior, 1 / (1 + (1 / 5) * ((1 - 5 / 38)^(-59) - 1)))
})

test_that("thompson2013 reproduces the two-year example (Section 2.3)", {
  d <- sighting_data(data.frame(year = 1:2, sightings = c(1, 0)))
  run <- function(p, px, avg = "annealed") {
    thompson2013(d, test_year = 2, p = p, prior_extant = px,
                 prior_extinct = 1 - px, average = avg)
  }
  # P(X2) = 0.5: PQ 0.411, PA 0.407, PD 0.375
  expect_equal(round(run(list(c(0.6, 0.8)), 0.5, "quenched")$posterior, 3), 0.411)
  expect_equal(round(run(list(c(0.6, 0.8)), 0.5)$posterior, 3), 0.408)  # paper 0.407 (truncated)
  expect_equal(run(list(c(0.6, 0.8)), 0.5)$posterior, 0.62 / 1.52)
  expect_equal(round(run(0.6, 0.5)$posterior, 3), 0.375)
  # P(X2) = 0.75: PQ 0.676, PA 0.674 (PD printed as 0.743, a typo for 0.643)
  expect_equal(round(run(list(c(0.6, 0.8)), 0.75, "quenched")$posterior, 3), 0.676)
  expect_equal(round(run(list(c(0.6, 0.8)), 0.75)$posterior, 3), 0.674)
  expect_equal(round(run(0.6, 0.75)$posterior, 3), 0.643)
  # eq 13: B = <(1-p)p> / <1-p>
  expect_equal(run(list(c(0.6, 0.8)), 0.5)$bayes_factor, (0.62 / 3) / 0.3)
  expect_true(is.na(run(list(c(0.6, 0.8)), 0.5, "quenched")$bayes_factor))
})

test_that("thompson2013 monk seal averages agree with the simulated values", {
  # eqs 21-25 are Monte Carlo estimates; the package integrates exactly.
  q <- function(r, ...) th_seal(p = list(r), average = "quenched", ...)$posterior
  a <- function(r, ...) th_seal(p = list(r), ...)$posterior
  expect_equal(round(q(c(0.86, 0.92)), 4), 0.0065)  # paper 0.0064
  expect_equal(round(a(c(0.86, 0.92)), 4), 0.0064)  # paper 0.0065
  expect_lt(abs(q(c(0.86, 1)) - 0.0949), 0.002)     # 0.0933
  expect_lt(abs(a(c(0.86, 1)) - 0.0246), 0.0011)    # 0.0236
  expect_equal(round(q(c(0, 0.92)), 4), 0.0004)     # paper 0.0005
  expect_equal(round(a(c(0, 0.92)), 4), 0.0036)     # paper 0.0034
  expect_equal(round(q(c(0.86, 0.92), prior_extant = 0.5, prior_extinct = 0.5), 4),
               0.0002)
  expect_equal(round(a(c(0.86, 0.92), prior_extant = 0.5, prior_extinct = 0.5), 4),
               0.0002)
  # p in (0, 1), eq 22: paper 0.0102 / 0.0179 by simulation, heavy-tailed;
  # exact values 0.0131 / 0.0150
  expect_equal(round(q(c(0, 1)), 4), 0.0131)
  expect_equal(round(a(c(0, 1)), 4), 0.0150)
})

test_that("thompson2013 reproduces the three-type example (Table 1, eqs 28-31)", {
  s <- rbind(c(1, 1, 0), c(0, 0, 1), c(1, 0, 1), c(0, 1, 0), c(0, 0, 0),
             c(0, 1, 0), c(0, 0, 1))
  recs <- lapply(1:3, function(a) sighting_data(data.frame(year = 1:7, n = s[, a])))
  p <- list(c(0.2, 0.4), c(0.4, 0.6), c(0.3, 0.6))
  q <- list(c(0.2, 0.7), c(0.1, 0.5))
  run <- function(avg, ...) {
    thompson2013(recs[[1]], recs[2:3], 7, p = p, q = q, average = avg, ...)$posterior
  }
  # eq 29 (P(X_T) = P_E = 0.5): PQ 0.0231, PA 0.0151 (simulated)
  expect_lt(abs(run("quenched", prior_extant = 0.5, prior_extinct = 0.5) - 0.0231), 3e-4)
  expect_lt(abs(run("annealed", prior_extant = 0.5, prior_extinct = 0.5) - 0.0151), 3e-4)
  # eq 30 (defaults T_N/T, 1/T): PQ 0.0631, PA 0.0442 (simulated)
  expect_lt(abs(run("quenched") - 0.0631), 3e-4)
  expect_lt(abs(run("annealed") - 0.0442), 1e-3)
  # eq 31: deterministic at the midpoints, 0.0373
  d <- thompson2013(recs[[1]], recs[2:3], 7, p = c(0.3, 0.5, 0.45),
                    q = c(0.45, 0.3))
  expect_equal(round(d$posterior, 4), 0.0373)
  # quadrature has converged
  expect_equal(run("quenched", nodes = 16), run("quenched"), tolerance = 1e-6)
})

test_that("thompson2013 handles impossible false sightings", {
  d <- sighting_data(data.frame(year = 1:3, n = c(1, 0, 0)))
  u <- sighting_data(data.frame(year = 1:3, n = c(0, 0, 1)))
  res <- thompson2013(d, u, 3, p = c(0.5, 0.5), q = 1)
  expect_equal(res$posterior, 1)
})

test_that("thompson2013 validates its inputs", {
  expect_error(th_seal(p = 1.2), "probability")
  expect_error(th_seal(p = list(c(0.9, 0.8))), "range")
  expect_error(th_seal(p = c(0.5, 0.5)), "one element")
  expect_error(th_seal(), "`p`")
  expect_error(th_seal(p = 0.9, prior_extant = 1), "prior_extant")
  expect_error(thompson2013(seal13(), test_year = 1952, p = 0.9), "later than")
  u <- sighting_data(data.frame(year = 2015, n = 1))
  expect_error(thompson2013(seal13(), u, 2011, p = c(0.9, 0.9), q = 0.9), "precede")
  frac <- sighting_data(data.frame(year = c(1.5, 3), n = 1))
  expect_error(thompson2013(frac, test_year = 5, p = 0.9), "whole")
})
