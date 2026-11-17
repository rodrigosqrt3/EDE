alaska <- function(taxon, drop_oldest = FALSE) {
  d <- alaska_megafauna[alaska_megafauna$taxon == taxon, ]
  if (drop_oldest) d <- d[-which.max(d$age_ka_bp), ]
  list(sd = sighting_data(data.frame(time = -d$age_ka_bp, count = 1)),
       error = d$sd_ka[order(-d$age_ka_bp)])
}

test_that("the dataset is Table 1 of Solow, Roberts & Robbirt (2006)", {
  expect_equal(dim(alaska_megafauna), c(50, 3))
  expect_equal(table(alaska_megafauna$taxon)[["mammoth"]], 25)
  expect_equal(min(alaska_megafauna$age_ka_bp[alaska_megafauna$taxon ==
                                                "mammoth"]), 11.50)
  expect_equal(max(alaska_megafauna$age_ka_bp[alaska_megafauna$taxon ==
                                                "horse"]), 17.65)
})

test_that("maximum likelihood estimates reproduce the paper", {
  # Paper: mammoths 11.51 / 14.09 (m = 25), horses 12.36 / 16.68 (n = 24), in
  # ka BP. The mammoth profile likelihood is flat (0.001 in log between 11.504
  # and 11.51), so 0.01 ka is the resolution of the published figure.
  m <- alaska("mammoth")
  fm <- solow2006(m$sd, m$error)
  expect_lt(abs(-fm$estimate - 11.51), 0.01)
  expect_lt(abs(-fm$gamma - 14.09), 0.01)
  h <- alaska("horse", drop_oldest = TRUE)
  fh <- solow2006(h$sd, h$error)
  expect_lt(abs(-fh$estimate - 12.36), 0.01)
  expect_lt(abs(-fh$gamma - 16.68), 0.01)
  expect_equal(fm$k, 25)
  expect_equal(fh$k, 24)
  expect_s3_class(fm, "ede_estimate")
  expect_equal(fm$interval_type, "two-sided")
})

test_that("the settlement example of Solow (1997) is reproduced", {
  # Paper: maximum of the profile log-likelihood -113.26 at 651 yr BP, other
  # end-point 553 yr BP, 0.95 interval 617-687 yr BP. The earliest end is the
  # one wanted, so the age enters as it is.
  fit <- solow2006(data.frame(time = nz_shells$age_bp,
                              error = nz_shells$sd_yr))
  expect_equal(nrow(nz_shells), 22)
  expect_equal(round(fit$estimate), 651)
  expect_equal(round(fit$gamma), 553)
  expect_equal(round(fit$loglik, 2), -113.26)
  expect_equal(round(fit$lower), 617)
  expect_equal(round(fit$upper), 687)
  # the interval is the set where the profile log-likelihood exceeds
  # -113.26 - 3.84 / 2 = -115.18
  expect_equal(fit$loglik - qchisq(0.95, 1) / 2, -115.18, tolerance = 1e-4)
})

test_that("a data frame and a sighting_data object give the same fit", {
  m <- alaska("mammoth")
  a <- solow2006(m$sd, m$error)
  d <- alaska_megafauna[alaska_megafauna$taxon == "mammoth", ]
  b <- solow2006(data.frame(time = -d$age_ka_bp, error = d$sd_ka))
  expect_equal(a$estimate, b$estimate, tolerance = 1e-8)
  expect_equal(a$upper, b$upper, tolerance = 1e-6)
  # fewer than three ages cannot be used
  expect_error(solow2006(data.frame(time = c(1, 2), error = 0.1)),
               "between 3")
  expect_error(solow2006(data.frame(time = 1:5)), "columns `time` and `error`")
  expect_error(solow2006(data.frame(time = 1:5, error = c(1, 1, 1, 1, 0))),
               "positive")
})

test_that("the estimate maximises the likelihood and the limits solve the LR equation", {
  m <- alaska("mammoth")
  fit <- solow2006(m$sd, m$error, alpha = 0.1)
  d <- fit$data
  ll <- function(th, ga) EDE:::sr_ll(th, ga, d$time, d$error)
  # no nearby point does better
  for (dt in c(-0.05, 0.05)) {
    expect_lte(ll(fit$estimate + dt, fit$gamma), fit$loglik + 1e-8)
  }
  # at the limits, the profile statistic equals the chi-square(1) quantile
  crit <- stats::qchisq(0.9, 1)
  full <- list(theta = fit$estimate, gamma = fit$gamma, ll = fit$loglik)
  for (lim in c(fit$lower, fit$upper)) {
    expect_equal(EDE:::sr_stat(lim, full, d$time, d$error), crit,
                 tolerance = 1e-5)
  }
  expect_lt(fit$lower, fit$estimate)
  expect_gt(fit$upper, fit$estimate)
})

test_that("a common shift and a count above one behave as expected", {
  m <- alaska("mammoth")
  a <- solow2006(m$sd, m$error)
  shifted <- m$sd
  shifted$time <- shifted$time + 100
  b <- solow2006(shifted, m$error)
  expect_equal(b$estimate, a$estimate + 100, tolerance = 1e-5)
  expect_equal(b$upper, a$upper + 100, tolerance = 1e-4)
  # k = all equals the default; fewer ages use the most recent ones
  c5 <- solow2006(m$sd, m$error, k = 10)
  expect_equal(c5$k, 10)
  expect_equal(c5$data$time, sort(-alaska_megafauna$age_ka_bp[
    alaska_megafauna$taxon == "mammoth"], decreasing = TRUE)[1:10])
})

test_that("the joint region reproduces the qualitative findings of the paper", {
  m <- alaska("mammoth")
  h <- alaska("horse", drop_oldest = TRUE)
  reg <- solow2006_region(m$sd, h$sd, m$error, h$error, n = 40)
  expect_s3_class(reg, "solow2006_region")
  g <- reg$grid[reg$grid$inside, ]
  # the maximum likelihood point is inside the region
  expect_true(any(g$statistic < 0.5))
  # paper: no value with beta1 > beta2 (mammoth extinct earlier than horse in
  # ka BP means theta_mammoth > theta_horse on the forward axis), but values
  # with the horse extinct after 12.0 ka BP (theta_horse > -12.0)
  expect_true(all(g$theta1 >= g$theta2))
  expect_gt(max(g$theta2), -12.0)
  expect_equal(reg$critical, stats::qchisq(0.95, 2))
  expect_output(print(reg), "Joint confidence region")
})

test_that("the simulation check returns a band around the 45 degree line", {
  m <- alaska("mammoth")
  fit <- solow2006(m$sd, m$error)
  set.seed(1)
  chk <- solow2006_check(fit, nsim = 200)
  expect_s3_class(chk, "solow2006_check")
  expect_equal(nrow(chk), 25)
  expect_true(all(chk$lower <= chk$upper))
  # the model fits these data: the observed ages sit inside the band
  expect_gte(mean(chk$observed >= chk$lower & chk$observed <= chk$upper), 0.9)
  expect_error(solow2006_check(1), "result of solow2006")
})

test_that("input is validated", {
  m <- alaska("mammoth")
  expect_error(solow2006(m$sd, 0), "positive")
  expect_error(solow2006(m$sd, c(0.1, 0.2)), "one per row")
  expect_error(solow2006(m$sd, m$error, k = 2), "between 3")
  expect_error(solow2006(m$sd, m$error, alpha = 2), "alpha")
  expect_error(solow2006(m$sd, m$error, k = 100), "between 3")
})

test_that("non-finite times and constant ages are refused", {
  expect_error(solow2006(data.frame(time = c(1, NA, 3), error = 0.1)),
               "finite numbers")
  expect_error(solow2006(data.frame(time = c(5, 5, 5), error = 0.1)),
               "must not all be equal")
})

test_that("a single error is recycled over the records", {
  m <- alaska("mammoth")
  a <- solow2006(m$sd, error = 0.12)
  b <- solow2006(m$sd, error = rep(0.12, nrow(m$sd)))
  expect_equal(a$estimate, b$estimate)
})

test_that("region and check validate their arguments and can be plotted", {
  m <- alaska("mammoth")
  h <- alaska("horse", drop_oldest = TRUE)
  fit <- solow2006(m$sd, m$error)
  expect_error(solow2006_region(m$sd, h$sd, m$error, h$error, n = 5),
               "at least 10")
  expect_error(solow2006_check(fit, nsim = 2), "at least 10")
  expect_error(solow2006_check(fit, level = 1), "level")
  reg <- solow2006_region(m$sd, h$sd, m$error, h$error, n = 12)
  set.seed(1)
  chk <- solow2006_check(fit, nsim = 20)
  expect_invisible(plot(reg))
  expect_invisible(plot(chk))
})

test_that("limits that do not exist are Inf, and the fit warns about it", {
  m <- alaska("mammoth")
  d <- data.frame(time = m$sd$time, error = m$error)
  # a profile statistic that never reaches the critical value
  testthat::local_mocked_bindings(sr_stat = function(...) 0)
  lim <- EDE:::sr_limits(list(theta = -11.5), d$time, d$error, 3.84)
  expect_equal(lim, c(-Inf, Inf))
  # a fit whose upper limit does not exist
  testthat::local_mocked_bindings(sr_limits = function(...) c(-12, Inf))
  expect_warning(res <- solow2006(d), "upper confidence limit does not exist")
  expect_identical(res$upper, Inf)
  # the joint region then extends each axis beyond the estimate
  testthat::local_mocked_bindings(
    sr_limits = function(fit, x, s, crit) c(fit$theta - 1, Inf)
  )
  h <- alaska("horse")
  reg <- solow2006_region(m$sd, h$sd, m$error, h$error, n = 10)
  expect_identical(reg$fit1$upper, Inf)
  expect_true(all(is.finite(reg$grid$theta1)))
  expect_gt(max(reg$grid$theta1), reg$fit1$estimate)
  expect_gt(max(reg$grid$theta2), reg$fit2$estimate)
})
