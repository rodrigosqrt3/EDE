# Builds the example sighting records bundled with EDE. Each record is
# transcribed from the published source cited in R/data.R.
# Run from the package root: source("data-raw/datasets.R")

# Roberts & Solow (2003), Nature 426: 245. Reliable sightings; 1628 has two.
dodo <- data.frame(
  year = c(1598, 1601, 1602, 1607, 1611, 1628, 1631, 1638, 1662),
  count = c(1, 1, 1, 1, 1, 2, 1, 1, 1)
)

# Solow (1993a), Ecology 74: 962-964.
monk_seal <- data.frame(
  year = c(1915, 1922, 1932, 1948, 1952),
  count = 1
)

# Solow (1993b), J. Math. Biol. 32: 79-82, Table 2; also Jaric & Ebenhard
# (2010), Table 2. One row per sighting month.
ferret <- data.frame(
  year = c(rep(1972, 4), rep(1973, 6), rep(1974, 2), rep(1975, 3),
           rep(1976, 3), 1977, 1979, 1981, 1981, rep(1982, 3), 1983,
           1984, 1984),
  month = c(6, 7, 8, 10, 5:10, 6, 7, 5, 8, 10, 5, 9, 10, 6, 6, 9, 10,
            2, 3, 7, 7, 7, 9),
  count = 1
)

# Solow & Roberts (2003), Ecology 84: 1329-1332. Collections of
# Angraecum aff. rutenbergianum in Mauritius.
orchid <- data.frame(
  year = c(1933, 1941, 1952, 1962, 1965, 1966, 1967, 1971, 1974, 1976, 1978),
  count = 1
)

# Jaric & Roberts (2014), Biodivers. Conserv. 23: 2801-2815, Table 2, from
# Roberts, Elphick & Reed (2010).
bird_records <- function(species, scientific, physical, expert, controversial) {
  data.frame(
    species = species,
    scientific_name = scientific,
    year = c(physical, expert, controversial),
    evidence = rep(c("physical", "expert", "controversial"),
                   c(length(physical), length(expert), length(controversial)))
  )
}
extinct_birds <- rbind(
  bird_records(
    "Eskimo Curlew", "Numenius borealis",
    physical = c(1900:1906, 1908, 1911:1914, 1918, 1932, 1962, 1963),
    expert = c(1926, 1929, 1930, 1933, 1937, 1939, 1945, 1947, 1948, 1950,
               1959:1961),
    controversial = c(1956, 1964, 1968, 1970, 1972:1974, 1976, 1977,
                      1980:1983, 1985, 1987, 1990, 1992, 1996, 2002, 2006)
  ),
  bird_records(
    "Ivory-billed Woodpecker", "Campephilus principalis principalis",
    physical = c(1897:1902, 1904:1910, 1913, 1914, 1917, 1924, 1925, 1932,
                 1935, 1938, 1939),
    expert = c(1911, 1916, 1920, 1921, 1923, 1926, 1929:1931, 1933, 1934,
               1936, 1937, 1941:1944),
    controversial = c(1946, 1948:1952, 1955, 1958, 1959, 1962, 1966:1969,
                      1971:1974, 1976, 1981, 1982, 1985:1988, 1999,
                      2004:2006)
  ),
  bird_records(
    "Nukupuu", "Hemignathus lucidus",
    physical = c(1837, 1838, 1888, 1891:1897, 1899),
    expert = 1879,
    controversial = c(1960, 1961, 1965, 1967, 1968, 1973:1976, 1978, 1979,
                      1981, 1983, 1985:1996)
  ),
  bird_records(
    "Oahu Alauahio", "Paroreomyza maculata",
    physical = c(1837, 1888, 1891:1893, 1896, 1897, 1901:1903, 1950, 1968),
    expert = c(1936, 1937, 1939, 1940, 1946:1949),
    controversial = c(1952, 1955:1958, 1960, 1961, 1964, 1966, 1968, 1969,
                      1972:1978, 1984, 1985, 1989:1991, 1993, 1996, 1997,
                      2000:2002)
  )
)

# Solow, Roberts & Robbirt (2006), PNAS 103: 7351-7353, Table 1 (taken by the
# authors from Guthrie 2003, 2004): the 25 most recent radiocarbon ages of
# Alaskan mammoths and caballoid horses, in ka BP, with their standard
# deviations (ka). The oldest horse age (17.65) is the one the authors omit.
alaska_megafauna <- data.frame(
  taxon = rep(c("mammoth", "horse"), each = 25),
  age_ka_bp = c(
    11.50, 11.54, 11.91, 12.12, 12.19, 12.34, 12.43, 12.44, 12.48, 12.49,
    12.51, 12.58, 12.68, 12.88, 13.06, 13.23, 13.29, 13.34, 13.38, 13.41,
    13.44, 13.66, 13.69, 14.02, 14.09,
    12.48, 12.51, 12.56, 12.58, 12.71, 12.86, 13.25, 13.27, 13.94, 14.00,
    14.12, 14.26, 14.30, 14.62, 14.99, 15.37, 15.50, 15.57, 15.75, 15.92,
    16.13, 16.15, 16.48, 16.70, 17.65
  ),
  sd_ka = c(
    .16, .14, .13, .09, .13, .11, .18, .13, .08, .17,
    .15, .15, .14, .12, .15, .09, .14, .15, .09, .15,
    .09, .16, .19, .10, .16,
    .08, .13, .14, .14, .17, .14, .16, .15, .16, .16,
    .18, .16, .16, .17, .19, .19, .10, .19, .19, .19,
    .24, .21, .21, .22, .23
  )
)

# Solow (1997), Radiocarbon 39: 351-354, Table 1: approximate calibrated ages
# (yr BP) and error standard deviations of 22 shell samples from New Zealand,
# read by the author from Figure 3 of Higham & Hogg (1997).
nz_shells <- data.frame(
  sample = c("Wk-2282", "Wk-2768", "Wk-2770", "Wk-1083", "Wk-2502", "Wk-1515",
             "Wk-2515", "Wk-1167", "Wk-1737", "Wk-2751", "Wk-2410", "Wk-2411",
             "Wk-2412", "Wk-2362", "Wk-2508", "Wk-2752", "Wk-2856", "Wk-2632",
             "Wk-2857", "Wk-2440", "Wk-2441", "Wk-2851"),
  age_bp = c(670, 670, 620, 530, 565, 595, 610, 580, 610, 565, 580, 595, 635,
             580, 640, 565, 575, 575, 560, 625, 645, 570),
  sd_yr = c(30, 30, 40, 25, 55, 50, 50, 35, 40, 45, 50, 50, 30, 45, 25, 45,
            40, 40, 45, 30, 25, 40)
)

# Solow (1996), Paleobiology 22: 406-410, Table 1: number of finds and
# approximate location (m) of the lowest and highest find of the ten uppermost
# ammonite taxa of Seymour Island, from Macellari (1986); taxa 3 and 8 from
# Strauss & Sadler (1989). Taxa are numbered as in the table.
seymour_ammonites <- data.frame(
  taxon = 1:10,
  n = c(14, 24, 16, 8, 10, 11, 10, 8, 3, 4),
  lowest = c(480, 640, 800, 900, 900, 925, 955, 979, 990, 1100),
  highest = c(1165, 1165, 1147, 1145, 1135, 1120, 1130, 1091, 1120, 1145)
)

save(dodo, file = "data/dodo.rda", compress = "xz")
save(monk_seal, file = "data/monk_seal.rda", compress = "xz")
save(ferret, file = "data/ferret.rda", compress = "xz")
save(orchid, file = "data/orchid.rda", compress = "xz")
save(extinct_birds, file = "data/extinct_birds.rda", compress = "xz")
save(alaska_megafauna, file = "data/alaska_megafauna.rda", compress = "xz")
save(nz_shells, file = "data/nz_shells.rda", compress = "xz")
save(seymour_ammonites, file = "data/seymour_ammonites.rda", compress = "xz")
