# Tasmanian red fox carcasses, Table 1 of Caley & Barry (2014); observation
# from 2001 (y = 1,0,1,0,1,1,0,...,0 to 2012) or from 1999 (two leading zeros)
fox01 <- function() {
  sighting_data(data.frame(year = c(2001, 2003, 2005, 2006), count = 1))
}
fox99 <- function() {
  sighting_data(data.frame(year = c(1999, 2001, 2003, 2005, 2006),
                           count = c(0, 1, 1, 1, 1)))
}
near <- function(x, y, tol) expect_lt(max(abs(x - y)), tol)

test_that("constant model reproduces the published results (2001 start)", {
  res <- caley_barry2014(fox01(), test_year = 2012)
  expect_s3_class(res, "ede_bayes")
  # Results: P(extant as of 2013) = 6.2%; phi median 18.6% (2.3-50.7);
  # lambda median 51.3% (20.3-81.9)
  expect_equal(round(res$posterior, 3), 0.063)
  near(res$posterior, 0.062, 0.001)
  expect_equal(round(100 * unname(res$phi), 1), c(18.6, 2.3, 50.7))
  near(unname(res$lambda), c(0.513, 0.203, 0.819), 0.001)
  # Fig. 3B: the first bar (labelled 2007) is 0.795; it pools 2007 and 2008
  expect_equal(round(sum(res$extinction_time$probability[1:2]), 3), 0.795)
  expect_equal(res$extinction_time$time, 2007:2012)
  expect_equal(sum(res$extinction_time$probability) + res$posterior, 1)
})

test_that("constant model reproduces the published results (1999 start)", {
  res <- caley_barry2014(fox99(), test_year = 2012)
  # 13.9%; phi 14.0% (1.3-42.0); lambda 40.1% (15.3-70.4); Fig. 3B 0.67
  near(res$posterior, 0.139, 0.001)
  expect_equal(round(100 * unname(res$phi), 1), c(14.0, 1.3, 42.0))
  expect_equal(round(100 * unname(res$lambda), 1), c(40.1, 15.2, 70.5))
  near(unname(res$lambda[2]), 0.153, 0.001)
  near(unname(res$lambda[3]), 0.704, 0.001)
  expect_equal(round(sum(res$extinction_time$probability[1:2]), 2), 0.67)
})

test_that("constant model reproduces Fig. 3A projections", {
  d01 <- caley_barry2014(fox01(), test_year = 2015, data_out = TRUE)
  d99 <- caley_barry2014(fox99(), test_year = 2015, data_out = TRUE)
  expect_named(d01, c("time", "bayes_factor", "posterior"))
  expect_equal(d01$time, 2007:2015)
  # text: 1.6% (2001 start) and 4.7% (1999 start) in 2015, from a 1e6
  # iteration chain
  near(d01$posterior[d01$time == 2015], 0.016, 0.001)
  near(d99$posterior[d99$time == 2015], 0.047, 0.001)
  # values read from Fig. 3A (to about 0.01)
  near(d01$posterior[1:5], c(0.75, 0.49, 0.30, 0.17, 0.105), 0.011)
  near(d99$posterior[1:5], c(0.82, 0.62, 0.44, 0.30, 0.20), 0.011)
  expect_true(all(diff(d01$posterior) < 0))
  expect_equal(d01$posterior[d01$time == 2012],
               caley_barry2014(fox01(), test_year = 2012)$posterior)
})

test_that("published formulation equals the equations with lambda ~ Beta(a, b + 1)", {
  a <- caley_barry2014(fox01(), 2012)
  b <- caley_barry2014(fox01(), 2012, formulation = "equations",
                       prior_lambda = c(1.1, 2.1))
  expect_equal(a$posterior, b$posterior)
  expect_equal(a$phi, b$phi)
})

test_that("constant model matches brute-force integration of the equations", {
  y <- c(1, 0, 1, 0, 1, 1, rep(0, 6))
  g <- (seq_len(800) - 0.5) / 800
  W <- outer(g, rep(1, 800))
  L <- outer(rep(1, 800), g)
  prior <- dbeta(W, 1.1, 1.1) * dbeta(L, 1.1, 1.1)
  ext <- 0
  for (z in 7:12) ext <- ext + (1 - W)^(z - 1) * W * L^4 * (1 - L)^(z - 1 - 4)
  alive <- (1 - W)^12 * L^4 * (1 - L)^8
  brute <- sum(alive * prior) / sum((ext + alive) * prior)
  res <- caley_barry2014(fox01(), 2012, formulation = "equations")
  near(res$posterior, brute, 1e-4)
  expect_equal(res$prior, beta(1.1, 13.1) / beta(1.1, 1.1))
  odds <- function(p) p / (1 - p)
  expect_equal(res$bayes_factor, odds(res$posterior) / odds(res$prior))
})

test_that("non-constant model reproduces the published results", {
  skip_on_cran()
  # Results: r median -0.28 (-0.72 to -0.04); d 1.92 (0.48-4.3);
  # e0 -13.7 (-19.7 to -3.3); e1 5.5 (0.3-18.5); P(extant 2013) 26.9%.
  # Tolerances reflect Monte Carlo error over seeds with 1e5 iterations.
  set.seed(20140430)
  res <- caley_barry2014(fox01(), 2012, model = "nonconstant",
                         n_iter = 1.1e5, burn_in = 1e4)
  near(res$posterior, 0.269, 0.04)
  near(res$growth[["median"]], -0.28, 0.03)
  near(res$growth[["lower"]], -0.72, 0.05)
  near(res$delta[["median"]], 1.92, 0.15)
  near(res$eps0[["median"]], -13.7, 1)
  near(res$eps1[["median"]], 5.5, 1.5)
  expect_true(is.na(res$prior))
  expect_equal(nrow(res$samples), 1e5)
  expect_equal(sum(res$extinction_time$probability) + res$posterior, 1)
  # 1999 start: r -0.12 (-0.37 to 0.05); d 0.88 (0.20-2.57); 40.7%
  set.seed(19990101)
  res <- caley_barry2014(fox99(), 2012, model = "nonconstant",
                         n_iter = 1.1e5, burn_in = 1e4)
  near(res$posterior, 0.407, 0.04)
  near(res$growth[["median"]], -0.12, 0.03)
  near(res$delta[["median"]], 0.88, 0.1)
})

test_that("non-constant model returns a consistent short chain", {
  # a short run that is not skipped on CRAN, so the sampler is exercised by
  # every check; the published values are tested above with long chains
  set.seed(1)
  res <- caley_barry2014(fox01(), 2012, model = "nonconstant",
                         n_iter = 600, burn_in = 100)
  expect_s3_class(res, "ede_bayes")
  expect_true(is.na(res$prior) && is.na(res$bayes_factor))
  expect_gte(res$posterior, 0)
  expect_lte(res$posterior, 1)
  near(res$posterior + sum(res$extinction_time$probability), 1, 1e-12)
  expect_equal(res$extinction_time$time, 2007:2012)
  expect_equal(nrow(res$samples), 500L)
  expect_true(all(res$acceptance >= 0 & res$acceptance <= 1))
  for (par in c("growth", "delta", "eps0", "eps1")) {
    expect_named(res[[par]], c("median", "lower", "upper"))
    expect_true(res[[par]][["lower"]] <= res[[par]][["upper"]])
  }
  expect_true(all(res$samples$growth > -2.3 & res$samples$growth < 0.69))
  expect_true(all(res$samples$delta > 0.01 & res$samples$delta < 4.6))
  known <- res$samples$extinction_time[!is.na(res$samples$extinction_time)]
  expect_true(all(known >= 2007 & known <= 2012))

  set.seed(1)
  again <- caley_barry2014(fox01(), 2012, model = "nonconstant",
                           n_iter = 600, burn_in = 100)
  expect_identical(again$samples, res$samples)
})

test_that("caley_barry2014 validates its inputs", {
  expect_error(caley_barry2014(fox01(), 2012, model = "nonconstant",
                               prior_delta = c(-1, 4.6)), "prior_delta")
  expect_error(caley_barry2014(fox01(), 2012, model = "nonconstant",
                               prior_eps0 = 1), "prior_eps0")
  expect_error(caley_barry2014(fox01()), "test_year")
  expect_error(caley_barry2014(fox01(), 2006), "later than")
  expect_error(caley_barry2014(fox01(), 2012.5), "whole number")
  frac <- sighting_data(data.frame(year = c(2001, 2003.5), count = 1))
  expect_error(caley_barry2014(frac, 2010), "whole-number")
  expect_error(caley_barry2014(fox01(), 2012, prior_phi = c(1, -1)), "prior_phi")
  expect_error(caley_barry2014(fox01(), 2012, prior_lambda = 1), "prior_lambda")
  expect_error(caley_barry2014(fox01(), 2012, model = "nonconstant",
                               data_out = TRUE), "constant model")
  expect_error(caley_barry2014(fox01(), 2012, model = "nonconstant",
                               prior_growth = c(0.5, 0.69)), "starting values")
  expect_error(caley_barry2014(fox01(), 2012, model = "nonconstant",
                               n_iter = 10, burn_in = 10), "n_iter")
  none <- sighting_data(data.frame(year = 2001, count = 0))
  expect_error(caley_barry2014(none, 2012), "no sightings")
})

test_that("print.ede_bayes works for caley_barry2014", {
  out <- capture.output(print.ede_bayes(caley_barry2014(fox01(), 2012)))
  expect_match(out[1], "Caley & Barry")
  expect_match(out[3], "0\\.06263")
})
