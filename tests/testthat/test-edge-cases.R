# Input validation and numerical edge cases across methods

no_sightings <- function() {
  sighting_data(data.frame(year = c(1900, 1910), count = 0))
}
simple_record <- function() {
  sighting_data(data.frame(year = c(1900, 1905, 1910, 1915, 1920), count = 1))
}

test_that("methods reject a record with no sightings", {
  expect_error(solow1993_bayes(no_sightings(), test_year = 1950), "no sightings")
  expect_error(solow1993b(no_sightings(), test_year = 1950), "no sightings")
  expect_error(rout2009(no_sightings(), test_year = 1950), "no sightings")
  expect_error(burgman1995(no_sightings(), test_year = 1950), "no sightings")
  expect_error(brook2019(no_sightings(), reliability = c(1, 1)), "no sightings")
})

test_that("expand_record rejects a test year before the last sighting", {
  expect_error(burgman1995(simple_record(), test_year = 1912), "later than")
})

test_that("burgman p-values warn when cancellation destroys precision", {
  expect_warning(EDE:::burgman_equation3_pvalue(1600, 1280, 1), "equation 3 lost precision")
  expect_warning(EDE:::burgman_equation4_pvalue(200, 200, 1), "equation 4 lost precision")
})

test_that("brook2019 validates its arguments and failing methods", {
  sd <- simple_record()
  expect_error(brook2019(sd, reliability = rep(1, 5), min_sightings = 0), "min_sightings")
  expect_error(brook2019(sd, reliability = list(1, 1)), "one element per")
  expect_error(
    brook2019(sd, reliability = rep(1, 5), method = function(sd, ...) c(1, 2),
              n_resamples = 5),
    "no usable value"
  )
  res <- brook2019(sd, reliability = rep(1, 5), n_resamples = 5,
                   persist_year = 1925, seed = 1)
  expect_equal(res$persist_year, 1925)
  expect_equal(res$persistence, mean(res$values > 1925))
})

test_that("marshall1994_trend has a print method", {
  sd <- sighting_data(data.frame(time = c(0, 1, 3, 4, 8, 9, 15), count = 1))
  out <- capture.output(print(marshall1994_trend(sd)))
  expect_match(out[1], "transposition test")
  expect_match(out[2], "gaps, T =")
  expect_match(out[3], "two-sided p-value")
})

test_that("marshall1997 passes on errors other than a divergent tail", {
  sd <- sighting_data(data.frame(time = c(0, 2, 5, 10), count = 1))
  boom <- function(t) {
    if (any(t > 10)) stop("recovery not defined beyond 10")
    rep(1, length(t))
  }
  expect_error(marshall1997(sd, recovery = boom), "not defined beyond 10")
})

test_that("weiss_marshall1999_prior accepts one record and checks its input", {
  r1 <- sighting_data(data.frame(time = 1:5, count = c(1, 0, 1, 1, 0)))
  expect_equal(weiss_marshall1999_prior(r1), c(c = 2, d = 1))
  expect_error(weiss_marshall1999_prior(list(r1, data.frame(time = 1, count = 1))),
               "sighting_data")
  all_present <- sighting_data(data.frame(time = 1:3, count = 1))
  expect_warning(weiss_marshall1999_prior(all_present), "zero parameter")
})

test_that("weiss_marshall1999 quantile handles weights summing below prob", {
  post <- list(w = c(0.6, 0.39), u_lo = c(0, 1), u_hi = c(1, Inf), l = 10)
  expect_equal(EDE:::weiss_marshall1999_quantile(post, 0.995), Inf)
  post_finite <- list(w = c(0.6, 0.39), u_lo = c(0, 1), u_hi = c(1, 3), l = 10)
  expect_equal(EDE:::weiss_marshall1999_quantile(post_finite, 0.995), 3)
})

test_that("thompson2013 validates its inputs", {
  sd <- simple_record()
  expect_error(thompson2013(data.frame(year = 1900, count = 1), test_year = 1930, p = 0.5),
               "definite")
  expect_error(thompson2013(sd, uncertain = list(data.frame(year = 1900, count = 1)),
                            test_year = 1930, p = 0.5), "uncertain")
  expect_error(thompson2013(no_sightings(), test_year = 1930, p = 0.5), "no definite")
  expect_error(thompson2013(sd, test_year = 1930.5, p = 0.5), "whole number")
})

test_that("thompson2017 validates its inputs and prints Monte Carlo bounds", {
  rec <- data.frame(year = c(2000, 2005), pci = 1)
  sur <- data.frame(year = 2008, eps = 0.5, pi = 0.5, pr = 0.9)
  none <- c(eps = 0, pi = 0, pr = 0)
  expect_error(thompson2017(list(year = 2000), sur, none, 1), "data frame")
  expect_error(thompson2017(data.frame(year = 2000.5, pci = 1), sur, none, 1),
               "whole numbers")
  expect_error(thompson2017(rec, sur, list(eps = c(0, 0, 0), pi = 0, pr = 0), 1),
               "one value")
  expect_error(thompson2017(rec, sur, none, 1, end_year = 2010.5), "whole number")

  res <- thompson2017(rec, sur, none, 1, end_year = 2010, n_sim = 200)
  out <- capture.output(print(res))
  expect_true(any(grepl("Monte Carlo", out)))
})

test_that("solow2012_bayes validates its inputs and priors", {
  certain <- sighting_data(data.frame(year = c(1900, 1910, 1920), count = 1))
  uncertain <- sighting_data(data.frame(year = c(1930, 1940), count = 1))
  expect_error(solow2012_bayes(certain, data.frame(year = 1930, count = 1), 1950),
               "sighting_data object or NULL")
  expect_error(solow2012_bayes(no_sightings(), NULL, 1950), "no certain")
  expect_error(solow2012_bayes(certain, uncertain, 1950,
                               tau_prior = function(t) rep(-1, length(t))),
               "non-negative")
  expect_error(solow2012_bayes(certain, uncertain, 1950,
                               tau_prior = function(t) rep(0, length(t))),
               "integrates to zero")
  expect_error(solow2012_bayes(certain, uncertain, 1950, tau_prior = "exponential",
                               exp_quantile = c(1930, 1)), "prob must be")
  expect_error(solow2012_bayes(certain, uncertain, 1950, tau_prior = "exponential",
                               exp_quantile = c(1930, 1e-300)), "too extreme")
})
