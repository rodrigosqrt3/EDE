rec <- function(years) sighting_data(data.frame(year = years, sightings = 1))
ex <- rec(c(1900, 1904, 1909, 1911, 1913, 1918, 1920, 1925, 1928, 1930, 1934))

test_that("Fisher's g distribution matches simulation", {
  set.seed(11)
  for (k in c(4, 12, 30)) {
    e <- matrix(rexp(1e5 * k), 1e5, k)
    g <- apply(e, 1, max) / rowSums(e)
    for (x in quantile(g, c(0.1, 0.5, 0.9))) {
      expect_equal(EDE:::fisher_g_upper(x, k), mean(g > x), tolerance = 0.01)
    }
  }
  expect_equal(EDE:::fisher_g_upper(1, 5), 0)
  expect_equal(EDE:::fisher_g_upper(0.1, 5), 1)
})

test_that("Hasofer-Wang W equals 1 / ((k - 1)(k G - 1))", {
  x <- c(10, 9.1, 8.7, 7.2, 6.9, 5)
  k <- length(x)
  xb <- mean(x)
  w <- k * (xb - x[k])^2 / ((k - 1) * sum((x - xb)^2))
  d <- x[-k] - x[k]
  g <- sum(d^2) / sum(d)^2
  expect_equal(w, 1 / ((k - 1) * (k * g - 1)))
  r <- EDE:::endpoint_hasofer_wang(x, k, "finite", 1e4)
  expect_equal(r$stat, w)
})

test_that("result structure and printing", {
  r <- ede_endpoint_test(ex, nsim = 1e4)
  expect_s3_class(r, "ede_endpoint_test")
  expect_equal(nrow(r$table), 2L)
  expect_named(r$table, c("test", "k", "statistic", "p_value", "p_asymptotic"))
  expect_true(all(r$table$p_value >= 0 & r$table$p_value <= 1))
  expect_output(print(r), "Gumbel domain")
})

test_that("statistics are invariant to origin and unit of time", {
  a <- ede_endpoint_test(ex, k = 6, nsim = 1e4)$table
  b <- ede_endpoint_test(rec(3 * c(1900, 1904, 1909, 1911, 1913, 1918, 1920,
                                   1925, 1928, 1930, 1934) + 50),
                         k = 6, nsim = 1e4)$table
  expect_equal(a$statistic, b$statistic)
  expect_equal(a$p_value, b$p_value)
})

test_that("one-sided p-values are complementary and two-sided is twice the min", {
  f <- ede_endpoint_test(ex, k = 6, test = "neves", alternative = "finite")
  i <- ede_endpoint_test(ex, k = 6, test = "neves", alternative = "infinite")
  tw <- ede_endpoint_test(ex, k = 6, test = "neves",
                          alternative = "two.sided")
  expect_equal(f$table$p_value + i$table$p_value, 1)
  expect_equal(tw$table$p_value, min(1, 2 * min(f$table$p_value,
                                                  i$table$p_value)))
})

test_that("size is controlled under an exponential tail", {
  skip_on_cran()
  set.seed(5)
  rej <- replicate(300, {
    s <- sort(rexp(50))
    ede_endpoint_test(rec(s), k = 8, test = "neves")$table$p_value < 0.05
  })
  expect_lt(mean(rej), 0.10)
})

test_that("a uniform sighting rate is flagged as finite", {
  skip_on_cran()
  set.seed(6)
  rej <- replicate(200, {
    s <- sort(runif(60))
    ede_endpoint_test(rec(s), k = 10, nsim = 1e4)$table$p_value < 0.05
  })
  expect_gt(min(rowMeans(rej)), 0.3)
})

test_that("input validation", {
  expect_error(ede_endpoint_test(1), "sighting_data")
  expect_error(ede_endpoint_test(ex, k = 2), "at least 4")
  expect_error(ede_endpoint_test(ex, k = 11, test = "neves"), "cannot exceed")
  expect_error(ede_endpoint_test(rec(c(1, 2, 3)), nsim = 1e4), "at least 4")
  expect_error(ede_endpoint_test(ex, nsim = 10), "nsim")
  expect_error(ede_endpoint_test(ex, alternative = "bad"))
})

# Fraga Alves, Neves & Rosario (2016): the G* statistic

times_desc <- function(sd) sort(sd$time[sd$count > 0], decreasing = TRUE)

rec20 <- sighting_data(data.frame(
  year = c(1890, 1894, 1899, 1901, 1903, 1908, 1910, 1915, 1918, 1920, 1924,
           1927, 1931, 1933, 1936, 1940, 1943, 1947, 1950, 1952),
  sightings = 1
))

test_that("G* is the statistic of equation 17, built on the estimator of (6)", {
  t <- times_desc(rec20)
  k <- 5
  est <- fraga_alves_neves2014(rec20, k = k)$estimate
  g <- (est - t[k + 1]) / (t[k + 1] - t[2 * k + 1])
  expect_equal(EDE:::rosario_stat(t, k), log(2) * g - (log(k) + log(2) / 2))
  res <- ede_endpoint_test(rec20, k = k, test = "rosario", nsim = 2000)
  expect_equal(res$table$statistic, EDE:::rosario_stat(t, k))
  expect_equal(res$table$test, "Fraga Alves, Neves & Rosario (2016)")
  expect_equal(res$k, k)
})

test_that("G* is invariant to origin and unit of time", {
  t <- times_desc(rec20)
  expect_equal(EDE:::rosario_stat(t, 4), EDE:::rosario_stat(10 * t + 123, 4))
})

test_that("the asymptotic p-value is the Gumbel distribution function", {
  res <- ede_endpoint_test(rec20, k = 4, test = "rosario", nsim = 2000)
  g <- res$table$statistic
  expect_equal(res$table$p_asymptotic, exp(-exp(-g)))
  inf <- ede_endpoint_test(rec20, k = 4, test = "rosario",
                           alternative = "infinite", nsim = 2000)
  expect_equal(inf$table$p_asymptotic, 1 - exp(-exp(-g)))
  two <- ede_endpoint_test(rec20, k = 4, test = "rosario",
                           alternative = "two.sided", nsim = 2000)
  expect_equal(two$table$p_asymptotic, min(1, 2 * min(exp(-exp(-g)),
                                                      1 - exp(-exp(-g)))))
})

test_that("G* has a Gumbel limit under an exponential tail", {
  # Theorem 6: for large k the null of G* approaches the Gumbel distribution
  # (mean 0.5772, standard deviation pi / sqrt(6)).
  null <- EDE:::rosario_null(60, 4000)
  expect_lt(abs(mean(null) - 0.5772), 0.25)
  expect_lt(abs(stats::sd(null) - pi / sqrt(6)), 0.25)
  # and the Gumbel p-value is roughly calibrated there
  q <- stats::quantile(null, c(0.1, 0.5, 0.9), names = FALSE)
  expect_lt(max(abs(exp(-exp(-q)) - c(0.1, 0.5, 0.9))), 0.1)
})

test_that("the simulated null is the exponential-tail null of the statistic", {
  # draw sightings from a standard exponential tail: the exact p-value is
  # uniform, so about 5% of records are rejected at 0.05
  set.seed(1)
  k <- 4
  p <- vapply(1:400, function(i) {
    x <- sort(stats::rexp(2 * k + 1), decreasing = TRUE)
    sd <- sighting_data(data.frame(time = x, count = 1))
    ede_endpoint_test(sd, k = k, test = "rosario", nsim = 4000)$table$p_value
  }, numeric(1))
  expect_lt(abs(mean(p <= 0.05) - 0.05), 0.04)
  expect_lt(abs(mean(p <= 0.5) - 0.5), 0.1)
})

test_that("it detects a finite end-point and an unbounded tail", {
  set.seed(2)
  k <- 10
  reject <- function(draw, alternative) {
    p <- vapply(1:40, function(i) {
      sd <- sighting_data(data.frame(time = sort(draw(100), decreasing = TRUE),
                                     count = 1))
      ede_endpoint_test(sd, k = k, test = "rosario", alternative = alternative,
                        nsim = 3000)$table$p_value
    }, numeric(1))
    mean(p <= 0.05)
  }
  # uniform sighting times: finite end-point (Weibull domain)
  expect_gt(reject(stats::runif, "finite"), 0.8)
  # Pareto with index 1: unbounded upper tail (Frechet domain); power 0.85 in
  # a longer run
  expect_gt(reject(function(n) 1 / stats::runif(n), "infinite"), 0.6)
})

test_that("it needs 2k + 1 sighting times and the default is unchanged", {
  expect_error(ede_endpoint_test(rec20, k = 10, test = "rosario"), "k")
  def <- ede_endpoint_test(rec20, k = 5, nsim = 2000)
  expect_equal(nrow(def$table), 2)
  both <- ede_endpoint_test(rec20, k = 5, test = c("neves", "rosario"),
                            nsim = 2000)
  expect_equal(both$table$test[1], "Neves, Picek & Fraga Alves (2006)")
  expect_equal(nrow(both$table), 2)
  expect_output(print(both), "Rosario")
  # default k respects the 2k + 1 limit when the test is requested
  auto <- ede_endpoint_test(rec20, test = "rosario", nsim = 2000)
  expect_lte(2 * auto$k + 1, 20)
})

test_that("the rosario test refuses sightings that do not span more than one time", {
  # 30 events at one time: the 4 latest sightings and the threshold coincide
  flat <- sighting_data(data.frame(time = c(0, 1), count = c(1, 30)))
  expect_error(ede_endpoint_test(flat, k = 4, test = "neves", events = TRUE),
               "span more than one")
  expect_error(ede_endpoint_test(flat, k = 4, test = "hasofer_wang",
                                 events = TRUE, nsim = 2000),
               "span more than one")
  expect_error(ede_endpoint_test(flat, k = 4, test = "rosario", events = TRUE,
                                 nsim = 2000),
               "span more than one")
})

test_that("the rosario null restores the global random stream", {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
  EDE:::rosario_null(3, 777)
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  # and when it did exist, it comes back unchanged
  set.seed(42)
  before <- get(".Random.seed", envir = globalenv())
  EDE:::rosario_null(3, 778)
  expect_identical(get(".Random.seed", envir = globalenv()), before)
})
