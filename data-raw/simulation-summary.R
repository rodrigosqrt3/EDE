# Summarise the raw output of data-raw/simulation.R
#
#   Rscript data-raw/simulation-summary.R [outdir]
#
# Reads <outdir>/simulation_raw.rds and writes <outdir>/simulation_summary.csv.
# The summary is also copied to inst/extdata/ for the vignette.
#
# Definitions (E = true extinction time)
#   ok            the method returned an estimate or an upper limit
#   failure       1 - mean(ok)
#   covered       lower <= E <= upper, a missing limit counting as unbounded
#   coverage      mean(covered) among ok fits
#   coverage_all  fits that cover E, divided by all replicates (a failure is a miss)
#   share_finite  share of ok fits whose upper limit is finite
#   excess_upper  median(upper - E) over fits with a finite upper limit
#   bias_median, mae_median  median of (estimate - E) and of |estimate - E|
# Coverage is reported only for methods that return an interval. A missing
# upper limit from `fraga_alves_neves2014()` and `dekkers1989()` means the
# method has no interval, not an unbounded one.

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1L) args[1] else "data-raw/simulation-output"
E <- 100
no_interval <- c("fraga_alves_neves2014", "dekkers1989")

raw <- readRDS(file.path(outdir, "simulation_raw.rds"))
raw <- raw[, c("scenario", "n", "rep", "method", "estimate", "lower", "upper")]

summarise_cell <- function(d) {
  ok <- !is.na(d$estimate) | !is.na(d$upper)
  has_int <- !(d$method[1] %in% no_interval)
  covered <- ok & (is.na(d$lower) | d$lower <= E) & (is.na(d$upper) | d$upper >= E)
  fin <- ok & is.finite(d$upper)
  data.frame(
    scenario = d$scenario[1], n = d$n[1], method = d$method[1], reps = nrow(d),
    failure = mean(!ok),
    bias_median = if (any(!is.na(d$estimate))) median(d$estimate - E, na.rm = TRUE) else NA_real_,
    mae_median = if (any(!is.na(d$estimate))) median(abs(d$estimate - E), na.rm = TRUE) else NA_real_,
    coverage = if (has_int && any(ok)) mean(covered[ok]) else NA_real_,
    coverage_all = if (has_int) sum(covered) / nrow(d) else NA_real_,
    share_finite = if (has_int && any(ok)) sum(fin) / sum(ok) else NA_real_,
    excess_upper = if (any(fin)) median(d$upper[fin] - E) else NA_real_,
    stringsAsFactors = FALSE
  )
}

summ <- do.call(rbind, lapply(
  split(raw, list(raw$scenario, raw$n, raw$method), drop = TRUE), summarise_cell))
rownames(summ) <- NULL
write.csv(summ, file.path(outdir, "simulation_summary.csv"), row.names = FALSE)
dir.create("inst/extdata", showWarnings = FALSE, recursive = TRUE)
file.copy(file.path(outdir, "simulation_summary.csv"),
          "inst/extdata/simulation_summary.csv", overwrite = TRUE)
cat("Wrote", nrow(summ), "rows.\n")
