# Yearly 0/1 records from the authors' code (Sighting_records.R of
# github.com/SarithaKodikara/Inferring-Extinction-Year-using-a-Bayesian-Approach),
# written as sighting years. Observation ends in 2010 for every species.
k20_ibw_c <- c(1898:1902, 1904:1910, 1913, 1914, 1917, 1924, 1925, 1932, 1935,
               1938, 1939)
k20_ibw_u <- c(1911, 1916, 1920, 1921, 1923, 1926, 1929:1931, 1933, 1934, 1936,
               1937, 1941:1944, 1946, 1948:1952, 1955, 1958, 1959, 1962,
               1966:1969, 1971:1974, 1976, 1981, 1982, 1985:1988, 1999,
               2004:2006)
k20_species <- list(
  nukupuu = list(start = 1838,
                 c = c(1838, 1888, 1891:1897, 1899),
                 u = c(1879, 1960, 1961, 1965, 1967, 1968, 1973:1976, 1978,
                       1979, 1981, 1983, 1985:1996)),
  eskimo = list(start = 1901,
                c = c(1901:1906, 1908, 1911:1914, 1918, 1932, 1962, 1963),
                u = c(1926, 1929, 1930, 1933, 1937, 1939, 1945, 1947, 1948,
                      1950, 1956, 1959:1961, 1964, 1968, 1970, 1972:1974, 1976,
                      1977, 1980:1983, 1985, 1987, 1990, 1992, 1996, 2002,
                      2006)),
  akialoa = list(start = 1888,
                 c = c(1888, 1891:1893, 1895, 1898:1900, 1960),
                 u = c(1896, 1897, 1936, 1941, 1957, 1963, 1965, 1969, 1992,
                       1995)),
  alauahio = list(start = 1837,
                  c = c(1888, 1891:1893, 1896, 1897, 1901:1903, 1950, 1968),
                  u = c(1936, 1937, 1939, 1940, 1946:1949, 1952, 1955:1958,
                        1960, 1961, 1964, 1966, 1969, 1972:1978, 1984, 1985,
                        1989:1991, 1993, 1996, 1997, 2000:2002))
)
# certain record with a zero-count row at the first year of observation
k20_sd <- function(years, start = min(years)) {
  all <- sort(unique(c(start, years)))
  sighting_data(data.frame(year = all, count = as.numeric(all %in% years)))
}
k20_fit <- function(sp, model) {
  s <- k20_species[[sp]]
  kodikara2020(k20_sd(s$c, s$start), k20_sd(s$u), test_year = 2010,
               model = model)
}

test_that("IBW Models 1 and 2 reproduce Tables 2 and 3", {
  m1 <- kodikara2020(k20_sd(k20_ibw_c), test_year = 2010)
  # Table 2: 95% HDI 1940-1944, median 1940; p(tau_E <= 2010) = 1 (text)
  expect_equal(unname(m1$extinction_year), c(1940, 1940, 1944))
  expect_lt(m1$posterior, 1e-6)
  expect_s3_class(m1, "ede_bayes")
  expect_equal(m1$model, 1)
  # Fig. 2b: pc = 0.5 with 95% HDI between 0.3 and 0.6
  expect_equal(round(unname(m1$p_certain), 1), c(0.5, 0.3, 0.6))
  # Table 3: 95% HDI 1940-1945, median 1940
  m2 <- kodikara2020(k20_sd(k20_ibw_c), k20_sd(k20_ibw_u), test_year = 2010)
  expect_equal(m2$model, 2)
  expect_equal(unname(m2$extinction_year), c(1940, 1940, 1945))
  expect_equal(sum(m2$extinction_time$probability) + m2$posterior, 1)
})

test_that("Table 4 is reproduced (exact vs. JAGS, Monte Carlo error)", {
  # Table 4: median, 95% HDI for Models 1-3 (JAGS, 40000 draws). Exact
  # values here; the HDI upper bound of the heavy tail of tau_E carries
  # the most Monte Carlo error.
  tab <- list(
    nukupuu = rbind(c(1903, 1900, 1916), c(1903, 1900, 1914),
                    c(1903, 1900, 1924)),
    eskimo = rbind(c(1966, 1964, 1974), c(1967, 1964, 1984),
                   c(1969, 1964, 1991)),
    akialoa = rbind(c(1965, 1961, 1988), c(1968, 1961, 1994),
                    c(1971, 1961, 2046)),
    alauahio = rbind(c(1979, 1969, 2651), c(2058, 1969, 3952),
                     c(2079, 1969, 4077))
  )
  for (sp in names(tab)) {
    for (m in 1:3) {
      got <- unname(k20_fit(sp, m)$extinction_year)
      exp_ <- tab[[sp]][m, ]
      expect_lte(abs(got[1] - exp_[1]), 1)               # median
      expect_equal(got[2], exp_[2])                      # HDI lower
      # HDI upper: within 2 years when it is before 2010. Beyond 2010 the
      # posterior tail is P(tau_E > K) ~ 1 / K, whose 95% point has a Monte
      # Carlo relative error of about 4.4 / sqrt(n_eff), i.e. 7-14% for
      # n_eff of 1000-4000 draws; allow 15% of its distance from the start.
      s0 <- k20_species[[sp]]$start
      expect_lte(abs(got[3] - exp_[3]), max(2, 0.15 * (exp_[3] - s0)))
    }
  }
  # exact agreement where the posterior is concentrated before 2010
  expect_equal(unname(k20_fit("nukupuu", 1)$extinction_year), c(1903, 1900, 1916))
  expect_equal(unname(k20_fit("nukupuu", 2)$extinction_year), c(1903, 1900, 1914))
  expect_equal(unname(k20_fit("eskimo", 1)$extinction_year), c(1966, 1964, 1974))
  expect_equal(unname(k20_fit("eskimo", 3)$extinction_year), c(1969, 1964, 1991))
  expect_equal(unname(k20_fit("akialoa", 1)$extinction_year), c(1965, 1961, 1988))
  expect_equal(unname(k20_fit("akialoa", 2)$extinction_year), c(1968, 1961, 1994))
})

test_that("O'ahu 'Alauahio p(tau_E <= 2010) matches the text (0.76, 0.36)", {
  # Section 4.2: 0.76 under Model 1, 0.36 under Model 2 (JAGS)
  expect_lt(abs(1 - k20_fit("alauahio", 1)$posterior - 0.76), 0.01)
  expect_lt(abs(1 - k20_fit("alauahio", 2)$posterior - 0.36), 0.01)
})

test_that("IBW with all sightings treated as certain gives median near 2080", {
  # Section 3, "Treating uncertain sightings as certain": median 2080, with
  # the last sighting in 2007 as stated in the text (the authors' yearly
  # record ends in 2006; with it the median is 2039).
  all_c <- sort(c(k20_ibw_c, k20_ibw_u, 2007))
  fit <- kodikara2020(k20_sd(all_c), test_year = 2010)
  expect_lte(abs(fit$extinction_year[["median"]] - 2080), 2)
})

test_that("closed forms agree with direct numerical integration", {
  # Model 2 likelihood (eq. 9) integrated numerically over (pc, puv, pui)
  yc <- c(1, 1, 0, 1, 0, 0, 0, 0, 0, 0)
  yu <- c(0, 0, 1, 0, 1, 1, 0, 1, 0, 1)
  T <- length(yc)
  fit <- kodikara2020(k20_sd(which(yc == 1)), k20_sd(which(yu == 1)),
                      test_year = T, model = 2)
  marg <- function(k) {
    k1 <- min(k, T + 1) - 1
    a <- sum(yu[seq_len(k1)]); nc <- sum(yc); Nu <- sum(yu)
    b <- k1 - nc - a; cc <- Nu - a; d <- T - k1 - cc
    f_pu <- function(puv) vapply(puv, function(v) {
      stats::integrate(function(pui) {
        pu <- v * (1 - pui) + pui * (1 - v) + v * pui
        pu^a * (1 - pu)^b * pui^cc * (1 - pui)^d
      }, 0, 1, rel.tol = 1e-10)$value
    }, numeric(1))
    beta(nc + 1, k1 - nc + 1) *
      stats::integrate(f_pu, 0, 1, rel.tol = 1e-10)$value
  }
  ks <- seq(5, T + 1)
  w <- vapply(ks, marg, numeric(1)) *
    ifelse(ks <= T, 1 / (ks * (ks + 1)), 1 / (T + 1))
  w <- w / sum(w)
  expect_equal(c(fit$extinction_time$probability, fit$posterior), w,
               tolerance = 1e-6)
  # Model 3 (puv = pc)
  fit3 <- kodikara2020(k20_sd(which(yc == 1)), k20_sd(which(yu == 1)),
                       test_year = T, model = 3)
  marg3 <- function(k) {
    k1 <- min(k, T + 1) - 1
    a <- sum(yu[seq_len(k1)]); nc <- sum(yc); Nu <- sum(yu)
    b <- k1 - nc - a; cc <- Nu - a; d <- T - k1 - cc
    f <- function(pc) vapply(pc, function(p) {
      stats::integrate(function(pui) {
        pu <- p * (1 - pui) + pui * (1 - p) + p * pui
        p^nc * ((1 - p) * pu)^a * ((1 - p) * (1 - pu))^b * pui^cc *
          (1 - pui)^d
      }, 0, 1, rel.tol = 1e-10)$value
    }, numeric(1))
    stats::integrate(f, 0, 1, rel.tol = 1e-10)$value
  }
  w3 <- vapply(ks, marg3, numeric(1)) *
    ifelse(ks <= T, 1 / (ks * (ks + 1)), 1 / (T + 1))
  w3 <- w3 / sum(w3)
  expect_equal(c(fit3$extinction_time$probability, fit3$posterior), w3,
               tolerance = 1e-6)
  expect_true(all(is.na(fit3$p_certain)))
})

test_that("prior, Bayes factor and output structure", {
  fit <- kodikara2020(k20_sd(c(1, 2, 4)), test_year = 9)
  expect_equal(fit$prior, 1 / 10)
  expect_equal(fit$bayes_factor,
               (fit$posterior / (1 - fit$posterior)) / (1 / 9))
  expect_equal(fit$extinction_time$time, 5:9)
  expect_equal(fit$level, 0.95)
  expect_output(print(fit), "Kodikara et al. \\(2020\\) Model 1")
  narrow <- kodikara2020(k20_sd(c(1, 2, 4)), test_year = 9, level = 0.5)
  expect_lte(narrow$extinction_year[["upper"]] - narrow$extinction_year[["lower"]],
             fit$extinction_year[["upper"]] - fit$extinction_year[["lower"]])
})

test_that("uncertain sightings in certain years are censored", {
  c_sd <- k20_sd(c(1, 2, 4))
  a <- kodikara2020(c_sd, k20_sd(c(3, 6)), test_year = 9, model = 2)
  b <- kodikara2020(c_sd, k20_sd(c(2, 3, 6)), test_year = 9, model = 2)
  expect_equal(a$extinction_time, b$extinction_time)
})

test_that("input validation", {
  c_sd <- k20_sd(c(1, 2, 4))
  expect_error(kodikara2020(data.frame(t = 1, c = 1), test_year = 5),
               "sighting_data")
  expect_error(kodikara2020(c_sd, uncertain = 1, test_year = 5),
               "sighting_data")
  expect_error(kodikara2020(c_sd, test_year = 4), "later than")
  expect_error(kodikara2020(c_sd, test_year = 9, model = 2), "uncertain")
  expect_error(kodikara2020(c_sd, test_year = 9, model = 4), "1, 2 or 3")
  expect_error(kodikara2020(c_sd, test_year = 9, level = 1), "level")
  expect_error(kodikara2020(c_sd, test_year = 9.5), "whole")
  expect_error(kodikara2020(k20_sd(c(1.5, 2)), test_year = 9), "whole")
  expect_error(kodikara2020(sighting_data(data.frame(t = 1:3, c = 0)),
                            test_year = 9), "no certain")
  expect_error(kodikara2020(c_sd, k20_sd(c(0, 6)), test_year = 9),
               "between")
  expect_error(kodikara2020(c_sd, k20_sd(c(3, 10)), test_year = 9),
               "between")
})
