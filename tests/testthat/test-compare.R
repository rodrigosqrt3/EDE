test_that("ede_compare collects every method into one table", {
  sd <- sighting_data(monk_seal)
  out <- ede_compare(sd, test_year = 2050)
  expect_named(out, c("method", "function_name", "estimate", "lower", "upper",
                      "interval_type", "note"))
  expect_equal(nrow(out), 10L)
  expect_equal(out$upper[out$function_name == "solow1993"],
               solow1993(sd, test_year = 2050)$upper)
  expect_equal(out$estimate[out$function_name == "ole"], ole(sd)$estimate)
})

test_that("ede_compare keeps methods that fail or warn, with a note", {
  # two sightings: OLE and the Weibull test need three
  sd <- sighting_data(data.frame(year = c(1900, 1910), count = 1))
  expect_no_warning(ede_compare(sd, test_year = 1915))
  out <- ede_compare(sd, test_year = 1915)
  expect_true(is.na(out$estimate[out$function_name == "ole"]))
  expect_match(out$note[out$function_name == "ole"], "at least 3")
  # Solow (1993a) does not reach alpha by 1915
  expect_match(out$note[out$function_name == "solow1993"], "never falls")
})

test_that("ede_compare validates its inputs", {
  sd <- sighting_data(monk_seal)
  expect_error(ede_compare(sd), "test_year")
  expect_error(ede_compare(sd, test_year = 1950), "later than")
  expect_error(ede_compare(sd, test_year = 2000, alpha = 2), "must be in")
})

test_that("bundled datasets reproduce their published examples", {
  expect_equal(round(ole(sighting_data(dodo))$estimate), 1690)

  seal <- solow1993(sighting_data(monk_seal), test_year = 1992, data_out = TRUE)
  expect_equal(round(seal$chance[seal$time == 1992], 3), 0.053)

  months <- (ferret$year - 1972) * 12 + ferret$month
  fer <- sighting_data(data.frame(time = c(0, months), count = c(0, ferret$count)))
  p <- solow1993b(fer, test_year = 229, data_out = TRUE)
  expect_equal(round(p$chance[p$time == 229], 3), 0.050)

  orch <- solow_roberts2003(sighting_data(orchid), test_year = 1996, data_out = TRUE)
  expect_equal(orch$chance[orch$time == 1996], 0.10)

  expect_setequal(unique(extinct_birds$species),
                  c("Eskimo Curlew", "Ivory-billed Woodpecker", "Nukupuu",
                    "Oahu Alauahio"))
  expect_setequal(unique(extinct_birds$evidence),
                  c("physical", "expert", "controversial"))
})
