# Sighting records and reliability classes from the authors' R code
# (Brook et al. 2019, Data S1; Zenodo doi:10.5281/zenodo.3247189,
# ede_examples.R): physical, expert and controversial sightings.
brook_records <- list(
  NP = list(pe = c(1845, 1854, 1872, 1880, 1894, 1912),
            ie = c(1865, rep(1929, 4), 1935, 1943, 1948),
            cs = c(1897, 1901, 1909, 1914, 1954, 1959, 1963, rep(1969, 3), 1970,
                   1972, 1978, 1979, 1979, 1980, 1980, 1981, 1982, 1983)),
  SB = list(pe = c(1842, 1882, 1889), ie = c(1843, 1866, 1883),
            cs = c(1897, 1909)),
  BC = list(pe = c(1845, 1922, 1924, 1978, 1998, 2002, 2004), ie = 1995,
            cs = 2009),
  AG = list(pe = c(1929, 1960, 1963, 1969), ie = c(1972, 1982, 1985),
            cs = c(1947, 1970, 1971, 1986, 1988))
)

# Build a sighting_data object plus a list of per-sighting reliabilities.
brook_build <- function(rec, rel = c(0.95, 0.5, 0.1)) {
  yr <- c(rec$pe, rec$ie, rec$cs)
  pr <- c(rep(rel[1], length(rec$pe)), rep(rel[2], length(rec$ie)),
          rep(rel[3], length(rec$cs)))
  u <- sort(unique(yr))
  list(sd = sighting_data(data.frame(time = u, count = as.numeric(table(factor(yr, levels = u))))),
       rel = lapply(u, function(y) pr[yr == y]))
}

test_that("brook2019 reproduces Table 1(b) of Brook et al. (2019)", {
  skip_on_cran()
  # Table 1(b): BBJ extinction year (the upper 95% bound) and probability of
  # persistence, with the base reliabilities of the authors' code
  # (physical 0.95, expert 0.5, controversial 0.1). Stochastic: loose
  # tolerances of 4 years and 0.04.
  table1 <- data.frame(
    sp = c("NP", "NP", "SB", "SB", "BC", "BC", "AG", "AG"),
    model = rep(c("S93", "OLE"), 4),
    year = c(1998, 2029, 1925, 1936, 2035, 2014, 1999, 2002),
    persist_year = rep(c(1989, 1960, 2014, 2010), each = 2),
    pp = c(0.565, 0.648, 0.000, 0.013, 0.999, 0.075, 0.001, 0.017)
  )
  for (i in seq_len(nrow(table1))) {
    b <- brook_build(brook_records[[table1$sp[i]]])
    m <- if (table1$model[i] == "S93") strauss1989 else ole
    res <- suppressMessages(brook2019(b$sd, b$rel, method = m, n_resamples = 4000,
                                      persist_year = table1$persist_year[i],
                                      seed = 2019))
    info <- paste(table1$sp[i], table1$model[i])
    expect_lt(abs(floor(res$upper) - table1$year[i]), 4.5, label = paste(info, "year"))
    expect_lt(abs(res$persistence - table1$pp[i]), 0.04, label = paste(info, "persistence"))
  }
})

test_that("brook2019 with all reliabilities 1 gives the plain method", {
  d <- sighting_data(data.frame(year = c(1915, 1922, 1932, 1948, 1952), count = 1))
  res <- brook2019(d, rep(1, 5), method = strauss1989, n_resamples = 50, seed = 1)
  plain <- strauss1989(d)
  expect_true(all(res$values == plain$estimate))
  expect_equal(res$estimate, plain$estimate)
  expect_equal(res$upper, plain$estimate)
  expect_equal(res$n_failed, 0)
  up <- brook2019(d, rep(1, 5), method = strauss1989, statistic = "upper",
                  n_resamples = 20, seed = 1)
  expect_equal(up$estimate, plain$upper)
})

test_that("brook2019 combines sightings within a time unit", {
  d <- sighting_data(data.frame(year = c(1900, 1910, 1920, 1930), count = c(1, 2, 1, 1)))
  res <- suppressMessages(brook2019(d, list(1, c(0.5, 0.5), 1, 0.2),
                                    n_resamples = 10, seed = 1))
  expect_equal(res$combined_reliability, c(1, 0.75, 1, 0.2))
  expect_error(brook2019(d, list(1, 0.5, 1, 0.2)), "one value per sighting")
})

test_that("brook2019 keeps reliability-1 sightings and drops reliability-0 ones", {
  d <- sighting_data(data.frame(year = c(1900, 1910, 1920, 1930, 1940), count = 1))
  # the 1940 sighting is never kept, so LAD-like statistics stay below it
  lad <- function(sd) max(sd$time[sd$count > 0])
  res <- brook2019(d, c(1, 1, 1, 0.5, 0), method = lad, n_resamples = 200, seed = 3)
  expect_true(all(res$values %in% c(1920, 1930)))
  expect_lt(abs(mean(res$values == 1930) - 0.5), 0.15)
})

test_that("brook2019 keeps zero-count rows as the observation origin", {
  d <- sighting_data(data.frame(year = c(1900, 1910, 1920, 1930), count = c(0, 1, 1, 1)))
  origin <- function(sd) min(sd$time)
  res <- brook2019(d, c(0.5, 0.5, 1), method = origin, n_resamples = 50,
                   min_sightings = 1, seed = 1)
  expect_true(all(res$values == 1900))
})

test_that("brook2019 records and reports failed resamples", {
  d <- sighting_data(data.frame(year = c(1900, 1910, 1920, 1930), count = 1))
  expect_message(
    res <- brook2019(d, c(1, 1, 0.5, 0.5), n_resamples = 400, seed = 1),
    "discarded"
  )
  # fewer than 3 kept sightings with probability 0.25
  expect_lt(abs(res$prop_failed - 0.25), 0.07)
  expect_equal(res$n_failed + length(res$values), 400)

  failing <- function(sd) stop("nope")
  expect_error(suppressMessages(brook2019(d, rep(1, 4), method = failing, n_resamples = 5)),
               "no usable value")
})

test_that("brook2019 seed is reproducible and restores the RNG state", {
  d <- brook_build(brook_records$AG)
  set.seed(42)
  before <- .Random.seed
  a <- suppressMessages(brook2019(d$sd, d$rel, n_resamples = 100, seed = 7))
  expect_identical(.Random.seed, before)
  b <- suppressMessages(brook2019(d$sd, d$rel, n_resamples = 100, seed = 7))
  expect_identical(a$values, b$values)
  expect_s3_class(a, "ede_estimate")
  expect_match(a$method, "Brook et al. \\(2019\\).*Strauss")
})

test_that("brook2019 validates its inputs", {
  d <- sighting_data(data.frame(year = c(1900, 1910, 1920), count = 1))
  expect_error(brook2019(d), "reliability")
  expect_error(brook2019(d, c(1, 2, 1)), "\\[0, 1\\]")
  expect_error(brook2019(d, c(1, 1)), "one value per row")
  expect_error(brook2019(d, rep(1, 3), n_resamples = 0), "n_resamples")
  expect_error(brook2019(d, rep(1, 3), alpha = 1), "alpha")
  expect_error(brook2019(d, rep(1, 3), persist_year = "a"), "persist_year")
})
