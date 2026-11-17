dates_of <- function(taxon, drop_oldest = FALSE) {
  d <- alaska_megafauna[alaska_megafauna$taxon == taxon, ]
  if (drop_oldest) d <- d[-which.max(d$age_ka_bp), ]
  # years, increasing towards the present
  data.frame(time = -d$age_ka_bp * 1000, error = d$sd_ka * 1000)
}

test_that("the McInerny estimates for each k follow equation 3", {
  m <- dates_of("mammoth")
  res <- bradshaw2012(sighting_data(data.frame(time = m$time, count = 1)))
  tk <- res$theta_k
  # k = 2: records at 11,500 and 11,540 BP, a span of 40 years
  expect_equal(tk$k[1], 2)
  expect_equal(tk$span[1], 40)
  expect_equal(tk$theta[1], -11500 + log(0.05) / log(1 - 2 / 40))
  # k = 3: the third record is at 11,910 BP, a span of 410 years
  expect_equal(tk$theta[2], -11500 + log(0.05) / log(1 - 3 / 410))
  # Figure 2a of the paper: about 11,450 and 11,090 years BP
  expect_lt(abs(-tk$theta[1] - 11450), 20)
  expect_lt(abs(-tk$theta[2] - 11090), 20)
  expect_equal(nrow(tk), 24)
})

test_that("the inverse-weighted estimate is the weighted mean of equation 4", {
  m <- dates_of("mammoth")
  res <- bradshaw2012(sighting_data(data.frame(time = m$time, count = 1)))
  w <- 1 / res$theta_k$span
  expect_equal(res$theta_w, sum(w * res$theta_k$theta) / sum(w))
  expect_equal(res$estimate, res$theta_w)
  expect_true(is.na(res$lower) && is.na(res$upper))
  # Figure 2: the weighted estimates lie between the two bounds of the axes
  expect_gt(-res$theta_w, 11200)
  expect_lt(-res$theta_w, 11400)
  h <- dates_of("horse", TRUE)
  resh <- bradshaw2012(sighting_data(data.frame(time = h$time, count = 1)))
  expect_gt(-resh$theta_w, 12210)
  expect_lt(-resh$theta_w, 12415)
})

test_that("the GRIWM medians match the extinction dates read from the paper", {
  # The paper reports the extinctions of horses and mammoths centred on 12,200
  # and 11,150 14C years BP (Figure 4), with the 25 most recent ages.
  set.seed(2012)
  h <- bradshaw2012(dates_of("horse", TRUE), nsim = 3000, level = 0.95)
  m <- bradshaw2012(dates_of("mammoth"), nsim = 3000)
  expect_lt(abs(-h$estimate - 12200), 30)
  expect_lt(abs(-m$estimate - 11150), 30)
  expect_lt(h$lower, h$estimate)
  expect_gt(h$upper, h$estimate)
  expect_equal(h$interval_type, "two-sided")
  expect_equal(h$alpha, 0.05)
  # the horse series is estimated to end before the mammoth series, as claimed
  expect_lt(h$estimate, m$estimate)
  # the mammoth interval is wider than without dating error: no interval at all
  expect_s3_class(h, "ede_estimate")
  expect_match(h$method, "Bradshaw")
})

test_that("a data frame and a sighting_data object with errors agree", {
  m <- dates_of("mammoth")[1:8, ]
  a <- bradshaw2012(m, nsim = 10)
  expect_equal(a$theta_w,
               bradshaw2012(sighting_data(data.frame(time = m$time,
                                                     count = 1)))$theta_w)
  set.seed(1)
  s1 <- bradshaw2012(m, nsim = 50)$draws
  expect_length(s1, 50)
})

test_that("undefined estimates are left out and reported", {
  # yearly records: k records in a span of k years or fewer give no estimate
  sd <- sighting_data(data.frame(time = 1:6, count = 1))
  expect_error(bradshaw2012(sd), "undefined")
  # most recent first: 40, 39, 38, 0; spans 1 and 2 for k = 2 and 3 are not
  # larger than k, so only k = 4 has an estimate
  sd2 <- sighting_data(data.frame(time = c(0, 38, 39, 40), count = 1))
  res <- bradshaw2012(sd2)
  expect_equal(sum(is.na(res$theta_k$theta)), 2)
  expect_equal(res$theta_w, res$theta_k$theta[3])
  expect_true(is.finite(res$theta_w))
})

test_that("k limits the series to the most recent records", {
  m <- dates_of("mammoth")
  a <- bradshaw2012(m, nsim = 20, k = 10)
  expect_equal(nrow(a$theta_k), 9)
})

test_that("input is validated", {
  m <- dates_of("mammoth")
  expect_error(bradshaw2012(m, k = 1), "between 2")
  expect_error(bradshaw2012(m, k = 99), "between 2")
  expect_error(bradshaw2012(m, nsim = 0), "nsim")
  expect_error(bradshaw2012(m, level = 1), "level")
  expect_error(bradshaw2012(m, alpha = 0), "alpha")
  expect_error(bradshaw2012(data.frame(time = 1:5)), "error")
})

test_that("too few records and undefined resamples are reported", {
  expect_error(bradshaw2012(sighting_data(data.frame(time = 5, count = 1))),
               "at least 2")
  # two records one time unit apart: k = 2 never has a span above 2
  expect_error(bradshaw2012(data.frame(time = c(0, 1), error = 0.01),
                            nsim = 20),
               "undefined for every resampled")
})

test_that("resamples without an estimate are dropped and counted", {
  set.seed(1)
  expect_warning(
    res <- bradshaw2012(data.frame(time = c(0, 5), error = 3), nsim = 200),
    "gave no estimate"
  )
  expect_gt(res$failed, 0)
  expect_equal(length(res$draws) + res$failed, 200)
})
