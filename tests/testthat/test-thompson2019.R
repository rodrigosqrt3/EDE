# Alaotra Grebe from 1990, Table 1 of Thompson et al. (2019)
grebe19_surveys <- function() {
  y <- c(1990, 1993, 1994, 1997, 1998, 1999, 2000, 2004, 2009)
  data.frame(year = y,
             eps_lower = ifelse(y == 2000, 0.70, 0.80),
             eps_upper = ifelse(y == 2000, 0.90, 0.95),
             pi_lower = 0.90, pi_upper = 0.95, pr_lower = 0.70, pr_upper = 0.90)
}
grebe19_passive <- c(eps_lower = 0, eps_upper = 0.05, pi_lower = 0.10,
                     pi_upper = 0.65, pr_lower = 0.40, pr_upper = 0.60)

test_that("thompson2019 reproduces the Alaotra Grebe case study", {
  d <- thompson2019(NULL, grebe19_surveys(), grebe19_passive, p_initial = 0.5,
                    end_year = 1994, data_out = TRUE)
  # Section 4: b(u) = 2.837, b(u') = 1.0047 and cumulative B from 1990 to
  # 1994 of 2.837, 2.850, 2.864, 8.124, 23.05
  expect_equal(round(d$b[d$type == "survey"][1], 3), 2.837)
  expect_equal(round(d$b[d$type == "passive"][1], 4), 1.0047)
  expect_equal(round(d$bayes_factor_extinction, 3)[1:4],
               c(2.837, 2.850, 2.864, 8.124))
  expect_equal(round(d$bayes_factor_extinction[5], 2), 23.05)
  res <- thompson2019(NULL, grebe19_surveys(), grebe19_passive,
                      p_initial = 0.5, end_year = 1994)
  # P(E) = 0.958 by 1995
  expect_equal(round(1 - res$posterior, 3), 0.958)
  expect_s3_class(res, "ede_bayes")
  expect_equal(res$bayes_factor, 1 / res$bayes_factor_extinction)
  expect_equal(res$test_year, 1994)
})

test_that("thompson2019 follows equations 7, 8 and 11", {
  rec <- data.frame(year = c(1, 3), pci_lower = c(0.2, 0.6), pci_upper = c(0.4, 0.8))
  sv <- data.frame(year = 4, eps = 0.5, pi = 0.8, pr = 0.5)
  d <- thompson2019(rec, sv, c(eps = 0.1, pi = 0.5, pr = 0.5),
                    p_initial = c(0.2, 0.6), end_year = 5, data_out = TRUE)
  b <- c(1 - 0.3, 1 / (1 - 0.025), 1 - 0.7, 1 / (1 - 0.2), 1 / (1 - 0.025))
  expect_equal(d$b, b)
  expect_equal(d$bayes_factor_extinction, cumprod(b))
  expect_equal(d$posterior, 1 / (1 + cumprod(b) * (1 / 0.4 - 1)))
  # a certain record fixes the probability at one (section 2.1)
  d <- thompson2019(data.frame(year = 1, pci = 1), NULL,
                    c(eps = 0.1, pi = 0.5, pr = 0.5), 0.5, end_year = 3,
                    data_out = TRUE)
  expect_equal(d$posterior, c(1, 1, 1))
})

test_that("thompson2019 prints with the ede_bayes method", {
  res <- thompson2019(NULL, grebe19_surveys(), grebe19_passive, 0.5,
                      end_year = 1994)
  out <- capture.output(print.ede_bayes(res))
  expect_match(out[1], "Thompson et al. \\(2019\\)")
  expect_match(out[3], "0\\.04159")
})

test_that("thompson2019 validates its inputs", {
  sv <- grebe19_surveys()
  expect_error(thompson2019(NULL, sv, grebe19_passive, p_initial = 1), "\\(0, 1\\)")
  expect_error(thompson2019(NULL, sv, grebe19_passive, p_initial = 0), "\\(0, 1\\)")
  expect_error(thompson2019(NULL, sv), "passive")
  expect_error(
    thompson2019(data.frame(year = 1, pci = 1),
                 data.frame(year = 2, eps = 1, pi = 1, pr = 1),
                 c(eps = 0, pi = 0, pr = 0), 0.5),
    "cannot both"
  )
})
