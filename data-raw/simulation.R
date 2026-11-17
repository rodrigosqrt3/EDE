# Simulation study: coverage and accuracy of the EDE estimators
# ------------------------------------------------------------------
# Usage (from the package root, in a terminal):
#
#   Rscript data-raw/simulation.R [reps] [cores] [outdir]
#
#   reps   replicates per scenario            (default 1000)
#   cores  worker processes                   (default: detectCores() - 1)
#   outdir directory for the results          (default: data-raw/simulation-output)
#
# A quick smoke test:  Rscript data-raw/simulation.R 5 2 data-raw/sim-test
#
# Design
#   * True extinction time E = 100 time units. Sighting times are drawn on
#     (0, E) from the scenario density below, rounded up to whole time units
#     and aggregated into annual counts (the structure of a real record).
#   * Every method sees the record up to the (unknown) extinction and is asked
#     for an upper bound at level 1 - alpha = 0.95 with test_year = 3 * E.
#   * Reported per scenario and method (see data-raw/simulation-summary.R):
#     failure rate, median bias and absolute error, coverage of E, and the
#     share and excess of finite upper limits.
#
# Scenarios (density of sighting times on (0, E), r = E - t)
#   uniform   constant rate                      f ~ 1        (EVI = -1)
#   beta2     rate vanishing linearly at E        f ~ r        (EVI = -1/2)
#   beta4     rate vanishing as r^3               f ~ r^3      (EVI = -1/4)
#   eff_dec   sighting effort falling linearly to 20% of its start
#   eff_inc   sighting effort rising linearly from 20% to its full value
# For `uniform`, `eff_dec` and `eff_inc` the density is positive at E, i.e. the
# non-regular case of the Weibull domain (EVI = -1) in which no endpoint
# estimator is consistent; they probe robustness, not efficiency.
#
# Results are written scenario by scenario, so an interrupted run resumes
# from where it stopped.

args <- commandArgs(trailingOnly = TRUE)
reps   <- if (length(args) >= 1L) as.integer(args[1]) else 1000L
cores  <- if (length(args) >= 2L) as.integer(args[2]) else max(1L, parallel::detectCores() - 1L)
outdir <- if (length(args) >= 3L) args[3] else "data-raw/simulation-output"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

E          <- 100
test_year  <- 3 * E
alpha      <- 0.05
sizes      <- c(10L, 20L, 50L)
scenarios  <- c("uniform", "beta2", "beta4", "eff_dec", "eff_inc")
seed_base  <- 20261113L
time_limit <- 60    # seconds per method and replicate

# --- data generation ------------------------------------------------------
draw_times <- function(n, scenario, E) {
  if (scenario %in% c("uniform", "beta2", "beta4")) {
    a <- c(uniform = 1, beta2 = 2, beta4 = 4)[[scenario]]
    return(E * (1 - runif(n)^(1 / a)))
  }
  w <- switch(scenario,
    eff_dec = function(t) 1 - 0.8 * t / E,
    eff_inc = function(t) 0.2 + 0.8 * t / E
  )
  out <- numeric(0)
  while (length(out) < n) {            # rejection sampling, w <= 1
    t <- runif(n, 0, E)
    out <- c(out, t[runif(n) < w(t)])
  }
  out[seq_len(n)]
}

make_record <- function(n, scenario, E) {
  yr <- pmax(1, ceiling(draw_times(n, scenario, E)))
  tab <- table(yr)
  sighting_data(data.frame(time = as.numeric(names(tab)),
                           count = as.integer(tab)))
}

# --- estimators -----------------------------------------------------------
method_calls <- function(sd) {
  list(
    ole                   = function() ole(sd, alpha = alpha),
    robson1964            = function() robson1964(sd, alpha = alpha),
    strauss1989           = function() strauss1989(sd, alpha = alpha),
    solow1993             = function() solow1993(sd, alpha = alpha, test_year = test_year),
    solow1993b            = function() solow1993b(sd, alpha = alpha, test_year = test_year),
    solow2005             = function() solow2005(sd, alpha = alpha, test_year = test_year),
    mcinerny2006          = function() mcinerny2006(sd, alpha = alpha, test_year = test_year),
    burgman1995           = function() burgman1995(sd, alpha = alpha, test_year = test_year),
    solow_roberts2003     = function() solow_roberts2003(sd, alpha = alpha, test_year = test_year),
    jaric2010             = function() jaric2010(sd, alpha = alpha, test_year = test_year),
    hall_wang1999         = function() hall_wang1999(sd, alpha = alpha),
    fraga_alves_neves2014 = function() fraga_alves_neves2014(sd),
    dekkers1989           = function() dekkers1989(sd, alpha = alpha),
    smith1985             = function() smith1985(sd, alpha = alpha)
  )
}

run_method <- function(f) {
  res <- tryCatch(
    suppressWarnings({
      setTimeLimit(elapsed = time_limit, transient = TRUE)
      on.exit(setTimeLimit(elapsed = Inf), add = TRUE)
      f()
    }),
    error = function(e) NULL
  )
  if (is.null(res)) return(c(estimate = NA_real_, lower = NA_real_, upper = NA_real_))
  c(estimate = res$estimate, lower = res$lower, upper = res$upper)
}

one_task <- function(task) {
  set.seed(seed_base + task$id)
  sd <- make_record(task$n, task$scenario, E)
  calls <- method_calls(sd)
  rows <- lapply(names(calls), function(m) {
    r <- run_method(calls[[m]])
    data.frame(scenario = task$scenario, n = task$n, rep = task$rep, method = m,
               estimate = r[["estimate"]], lower = r[["lower"]], upper = r[["upper"]],
               n_years = nrow(sd), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# --- cluster ----------------------------------------------------------------
pkg_root <- normalizePath(".")
cl <- parallel::makeCluster(cores)
on.exit(parallel::stopCluster(cl), add = TRUE)
parallel::clusterExport(cl, c("E", "test_year", "alpha", "seed_base", "time_limit",
                              "draw_times", "make_record", "method_calls",
                              "run_method", "one_task", "pkg_root"))
invisible(parallel::clusterEvalQ(
  cl, { suppressMessages(devtools::load_all(pkg_root, quiet = TRUE)); NULL }))

for (sc in scenarios) for (n in sizes) {
  f <- file.path(outdir, sprintf("%s_n%d.rds", sc, n))
  if (file.exists(f)) { cat("skip", basename(f), "\n"); next }
  t0 <- Sys.time()
  sc_id <- match(sc, scenarios) * 1e6 + n * 1e3
  tasks <- lapply(seq_len(reps), function(i)
    list(id = sc_id + i, scenario = sc, n = n, rep = i))
  res <- do.call(rbind, parallel::parLapplyLB(cl, tasks, one_task))
  saveRDS(res, f)
  cat(sprintf("%-8s n=%-3d %d reps  %.1f min\n", sc, n, reps,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

# --- summary ----------------------------------------------------------------
files <- list.files(outdir, pattern = paste0("^(", paste(scenarios, collapse = "|"),
                                            ")_n[0-9]+\\.rds$"),
                   full.names = TRUE)
saveRDS(do.call(rbind, lapply(files, readRDS)),
        file.path(outdir, "simulation_raw.rds"))
source("data-raw/simulation-summary.R")   # reads outdir from the command line
cat("\nDone. Files in", normalizePath(outdir), "\n")
