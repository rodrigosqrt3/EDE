#' Sighting record of the dodo
#'
#' Reliable sightings of the dodo (*Raphus cucullatus*) in Mauritius, used by
#' Roberts & Solow (2003) to estimate its extinction date with [ole()]. Two
#' sightings were made in 1628. A later report from 1674 was excluded by the
#' authors as unreliable.
#'
#' @format A data frame with 9 rows: `year` and `count` (number of sightings).
#' @source Roberts, D. L., & Solow, A. R. (2003). Flightless birds: When did
#'   the dodo become extinct? Nature, 426(6964), 245.
#' @examples
#' ole(sighting_data(dodo))
"dodo"

#' Sighting record of the Caribbean monk seal
#'
#' Reliable sightings of the Caribbean monk seal (*Monachus tropicalis*) used
#' as the example in Solow (1993a), Burgman et al. (1995) and Solow (2005).
#' There is no natural start of the observation period, so the first sighting
#' defines the origin.
#'
#' @format A data frame with 5 rows: `year` and `count`.
#' @source Solow, A. R. (1993a). Inferring extinction from sighting data.
#'   Ecology, 74(3), 962-964.
#' @examples
#' solow1993(sighting_data(monk_seal), test_year = 2000)
#' solow1993_bayes(sighting_data(monk_seal), test_year = 1992)
"monk_seal"

#' Sighting record of the black-footed ferret
#'
#' Monthly sightings of the black-footed ferret (*Mustela nigripes*) in South
#' Dakota from 1972 to 1984, one row per month with a sighting, used by
#' Solow (1993b) and Jaric & Ebenhard (2010).
#'
#' The two papers index months differently. Solow (1993b) counts months from
#' January 1972 (= 1), starts observation at 0 and tests to `T = 229`;
#' Jaric & Ebenhard (2010) count from June 1972 (= 1) and test to
#' December 1990 (`T = 223`).
#'
#' @format A data frame with 28 rows: `year`, `month` (1-12) and `count`.
#' @source Solow, A. R. (1993b). Inferring extinction in a declining
#'   population. Journal of Mathematical Biology, 32(1), 79-82.
#'
#'   Jaric, I., & Ebenhard, T. (2010). A method for inferring extinction based
#'   on sighting records that change in frequency over time. Wildlife Biology,
#'   16(3), 267-275.
#' @examples
#' # Solow (1993b): months from January 1972, observation starting at 0
#' months <- (ferret$year - 1972) * 12 + ferret$month
#' sd <- sighting_data(data.frame(time = c(0, months), count = c(0, ferret$count)))
#' tail(solow1993b(sd, test_year = 229, data_out = TRUE), 1)
"ferret"

#' Collection record of a Mauritian orchid
#'
#' Collections of *Angraecum* aff. *rutenbergianum*, one of the nine Mauritian
#' orchids analysed by Solow & Roberts (2003). With 1996 as the end of the
#' observation period, their nonparametric test gives a p-value of 0.10.
#'
#' @format A data frame with 11 rows: `year` and `count`.
#' @source Solow, A. R., & Roberts, D. L. (2003). A nonparametric test for
#'   extinction based on a sighting record. Ecology, 84(5), 1329-1332.
#' @examples
#' solow_roberts2003(sighting_data(orchid), test_year = 2030)
"orchid"

#' Sighting records of four possibly extinct birds, by evidence quality
#'
#' Sighting years of the Eskimo Curlew, the Ivory-billed Woodpecker (U.S.),
#' the Nukupu'u and the O'ahu 'Alauahio, classified by Roberts et al. (2010)
#' into physical evidence, independent expert opinion and controversial
#' sightings, as tabulated by Jaric & Roberts (2014). Jaric & Roberts assign
#' reliabilities of 0.8-0.9, 0.6-0.8 and 0.1-0.4 to the three classes. A year
#' can appear under more than one class.
#'
#' @format A data frame with columns `species`, `scientific_name`, `year` and
#'   `evidence` (`"physical"`, `"expert"` or `"controversial"`).
#' @source Jaric, I., & Roberts, D. L. (2014). Accounting for observation
#'   reliability when inferring extinction based on sighting records.
#'   Biodiversity and Conservation, 23(11), 2801-2815, Table 2.
#'
#'   Roberts, D. L., Elphick, C. S., & Reed, J. M. (2010). Identifying
#'   anomalous reports of putatively extinct species and why it matters.
#'   Conservation Biology, 24(1), 189-196.
#' @examples
#' curlew <- extinct_birds[extinct_birds$species == "Eskimo Curlew", ]
#' # keep the most reliable class for each year
#' curlew <- curlew[order(curlew$year, match(curlew$evidence,
#'   c("physical", "expert", "controversial"))), ]
#' curlew <- curlew[!duplicated(curlew$year), ]
#' reliability <- c(physical = 0.85, expert = 0.7, controversial = 0.25)
#' jaric_roberts2014(
#'   sighting_data(data.frame(year = curlew$year, count = 1)),
#'   reliability = reliability[curlew$evidence]
#' )
"extinct_birds"
