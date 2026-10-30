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

save(dodo, file = "data/dodo.rda", compress = "xz")
save(monk_seal, file = "data/monk_seal.rda", compress = "xz")
save(ferret, file = "data/ferret.rda", compress = "xz")
save(orchid, file = "data/orchid.rda", compress = "xz")
save(extinct_birds, file = "data/extinct_birds.rda", compress = "xz")
