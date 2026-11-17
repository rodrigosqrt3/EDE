test_that("curves carry p_value and the deprecated chance alias", {
  d <- sighting_data(data.frame(year = c(1900, 1905, 1911, 1920, 1930),
                                sightings = 1))
  curves <- list(
    solow1993(d, test_year = 1960, data_out = TRUE),
    solow1993b(d, test_year = 1960, data_out = TRUE),
    solow2005(d, test_year = 1960, data_out = TRUE),
    mcinerny2006(d, test_year = 1960, data_out = TRUE),
    burgman1995(d, test_year = 1960, data_out = TRUE),
    solow_roberts2003(d, test_year = 1960, data_out = TRUE),
    jaric2010(d, test_year = 1960, data_out = TRUE),
    jaric_roberts2014(d, reliability = rep(1, 5), test_year = 1960,
                      data_out = TRUE)
  )
  for (cv in curves) {
    expect_true(all(c("time", "p_value", "chance") %in% names(cv)))
    expect_identical(cv$p_value, cv$chance)
  }
})

test_that("as.data.frame turns an estimate into one row", {
  d <- sighting_data(data.frame(year = c(1900, 1905, 1911, 1920, 1930),
                                sightings = 1))
  df <- as.data.frame(ole(d))
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 1L)
  expect_named(df, c("method", "estimate", "lower", "upper", "interval_type",
                     "alpha"))
  expect_equal(df$estimate, ole(d)$estimate)
  pt <- as.data.frame(solow1993(d, test_year = 2000))
  expect_true(is.na(pt$estimate))
  expect_equal(pt$interval_type, "one-sided")
})

test_that("ede_compare has its own class and includes the extreme-value methods", {
  res <- ede_compare(sighting_data(monk_seal), test_year = 2050)
  expect_s3_class(res, "ede_compare")
  expect_s3_class(res, "data.frame")
  expect_true(all(c("hall_wang1999", "fraga_alves_neves2014") %in%
                    res$function_name))
  expect_equal(nrow(res), 12L)
  fan <- res[res$function_name == "fraga_alves_neves2014", ]
  expect_false(is.na(fan$estimate))
})

test_that("plot methods draw without error and return their argument", {
  pdf(NULL)
  on.exit(dev.off())
  sd <- sighting_data(dodo)
  expect_identical(plot(sd), sd)
  est <- ole(sd)
  expect_identical(plot(est), est)
  hw <- suppressWarnings(hall_wang1999(sd, k = 5))
  expect_identical(plot(hw), hw)
  pt <- solow2005(sd, test_year = 2000)
  expect_identical(plot(pt), pt)
  cmp <- ede_compare(sd, test_year = 2000)
  expect_identical(plot(cmp), cmp)
  set.seed(1)
  sel <- ede_select_k(sd, nsim = 49)
  expect_identical(plot(sel), sel)
  # a lower limit without an upper one
  lo_only <- new_ede_estimate(10, lower = 5, method = "x", alpha = 0.05)
  expect_identical(plot(lo_only), lo_only)
})

test_that("plot errors clearly when there is nothing to draw", {
  empty <- new_ede_estimate(NA_real_, method = "x", alpha = 0.05)
  expect_error(plot(empty), "nothing to plot")
})

test_that("S3 print methods return their argument invisibly", {
  d <- sighting_data(data.frame(year = c(1900, 1905, 1911, 1920, 1930),
                                sightings = 1))
  for (obj in list(ole(d), solow1993_bayes(d, test_year = 1960))) {
    ret <- withVisible(print(obj))
    expect_false(ret$visible)
    expect_identical(ret$value, obj)
  }
})

test_that("plot.ede_compare refuses a table without results", {
  empty <- structure(
    data.frame(method = "m", function_name = "m", estimate = NA_real_,
               lower = NA_real_, upper = NA_real_, interval_type = NA_character_,
               note = ""),
    class = c("ede_compare", "data.frame")
  )
  expect_error(plot(empty), "nothing to plot")
})
