rec <- function() {
  sighting_data(data.frame(
    year = c(1880, 1900, 1904, 1909, 1911, 1913, 1918, 1920, 1925, 1928, 1930,
             1934),
    sightings = 1
  ))
}

test_that("forward_stop implements the rule of G'Sell et al. (2016)", {
  expect_equal(forward_stop(c(0.001, 0.002, 0.5, 0.9), 0.1), 2L)
  expect_equal(forward_stop(c(0.9, 0.001, 0.001), 0.1), 0L)
  expect_equal(forward_stop(c(0.001, 0.001, 0.001), 0.1), 3L)
  # the rule takes the largest j, even if the running mean exceeds alpha
  # in between
  p <- c(0.001, 0.5, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
         0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
         0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001)
  expect_equal(forward_stop(p, 0.1), length(p))
})

test_that("ede_select_k returns a well-formed table", {
  set.seed(1)
  sel <- ede_select_k(rec(), nsim = 99)
  expect_s3_class(sel, "ede_select_k")
  tab <- sel$table
  expect_named(tab, c("k", "statistic", "shape", "p_value", "estimable",
                      "rejected"))
  expect_equal(tab$k, 12:4)
  expect_true(all(tab$p_value > 0 & tab$p_value <= 1))
  expect_true(all(tab$statistic >= 1 / (tab$k - 1)))   # Greenwood lower bound
  # rejections form a leading block; the selected k is the first not rejected
  j <- sum(tab$rejected)
  expect_equal(tab$rejected, seq_len(nrow(tab)) <= j)
  expect_equal(sel$k, if (j < nrow(tab)) tab$k[j + 1L] else NA_integer_)
})

test_that("ede_select_k keeps the largest k when the model fits", {
  skip_on_cran()
  set.seed(11)
  gam <- cumsum(stats::rexp(14))^(1 / 2)
  d <- sighting_data(data.frame(time = -gam, count = 1))
  expect_equal(ede_select_k(d, nsim = 199)$k, 14L)
})

test_that("ede_select_k is reproducible under set.seed and prints", {
  set.seed(5)
  a <- ede_select_k(rec(), nsim = 99, kmax = 8)
  set.seed(5)
  b <- ede_select_k(rec(), nsim = 99, kmax = 8)
  expect_equal(a$table, b$table)
  expect_equal(a$table$k, 8:4)
  out <- capture.output(ret <- print(a))
  expect_match(out[1], "^<Ordered goodness-of-fit")
  expect_match(out[2], "selected k:")
  expect_identical(ret, a)
})

test_that("ede_select_k validates its inputs", {
  expect_error(ede_select_k(data.frame(1)), "sighting_data")
  expect_error(ede_select_k(rec(), alpha = 0), "must be in \\(0, 1\\)")
  expect_error(ede_select_k(rec(), nsim = 5), "at least 19")
  expect_error(ede_select_k(rec(), kmin = 3), "at least 4")
  expect_error(ede_select_k(rec(), kmax = 40), "cannot exceed")
  expect_error(ede_select_k(rec(), kmin = 8, kmax = 6), "at least 8")
  few <- sighting_data(data.frame(year = c(1900, 1910, 1920), sightings = 1))
  expect_error(ede_select_k(few), "at least 4")
})

test_that("a flat top is refused and rejection at every k is reported", {
  flat <- sighting_data(data.frame(time = c(1, 2), count = c(1, 20)))
  expect_error(ede_select_k(flat, events = TRUE, nsim = 19),
               "do not span more than one time")
  t <- c(6, 21, 27, 38, 40, 45, 47, 50, 93, 94, 96, 97, 98, 100)
  set.seed(1)
  expect_warning(
    sel <- ede_select_k(sighting_data(data.frame(time = t, count = 1)),
                        alpha = 0.999, nsim = 49),
    "rejected for every"
  )
  expect_true(is.na(sel$k))
  expect_output(print(sel), "selected k: none")
})
