# Second simulation study: choice of k and false extinction declarations
# ----------------------------------------------------------------------
#   Rscript data-raw/simulation-extra.R <part> [reps] [cores] [outdir]
#
#   part   "k"       sensitivity of the extreme-value estimators to k
#          "extant"  species still extant: how often is extinction declared?
#          "both"
#   reps   replicates per scenario            (default 1000 for k, 2000 for extant)
#   cores  worker processes                   (default: detectCores() - 1)
#   outdir results directory                  (default: data-raw/simulation-output)
#
# Smoke test:  Rscript data-raw/simulation-extra.R both 5 2 data-raw/sim-test
#
# Part "k"
#   Same generator and scenarios as data-raw/simulation.R (true extinction
#   time E = 100; n = 10, 20, 50 sighting events). For ole(), hall_wang1999(),
#   fraga_alves_neves2014() and dekkers1989() the number k of latest
#   sightings is fixed at 5, 10, 20, or all, or chosen by ede_select_k()
#   (nsim = 199, alpha = 0.1), capped at the number of sightings available.
#   Question: does choosing k repair the under-coverage of ole() and make the
#   Hall-Wang upper limit finite more often?
#
# Part "extant"
#   The species persists. Sighting times are drawn on (0, 100), rounded up to
#   whole years, and every interval method is asked at test_year = 101 for
#   its 95% upper bound. A false alarm is an upper bound at or below 101,
#   i.e. extinction declared although the species is extant. Nominal rate for
#   a valid one-sided bound: at most 0.05.
#   Scenarios: constant rate; rate falling by a factor exp(2) (about 7)
#   over the window; rate falling by exp(3) (about 20); n = 5, 10, 20, 50.
#
# Both parts write one file per cell and resume when interrupted.

args <- commandArgs(trailingOnly = TRUE)
part   <- if (length(args) >= 1L) args[1] else "both"
reps_k <- if (length(args) >= 2L) as.integer(args[2]) else 1000L
reps_e <- if (length(args) >= 2L) as.integer(args[2]) else 2000L
cores  <- if (length(args) >= 3L) as.integer(args[3]) else max(1L, parallel::detectCores() - 1L)
outdir <- if (length(args) >= 4L) args[4] else "data-raw/simulation-output"
stopifnot(part %in% c("k", "extant", "both"))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

E <- 100; alpha <- 0.05; time_limit <- 120
seed_base <- 20261114L

# --- generators -------------------------------------------------------------
draw_times <- function(n, scenario, E) {
  if (scenario %in% c("uniform", "beta2", "beta4")) {
    a <- c(uniform = 1, beta2 = 2, beta4 = 4)[[scenario]]
    return(E * (1 - runif(n)^(1 / a)))
  }
  w <- switch(scenario,
    eff_dec = function(t) 1 - 0.8 * t / E,
    eff_inc = function(t) 0.2 + 0.8 * t / E,
    decl2   = function(t) exp(-2 * t / E),
    decl3   = function(t) exp(-3 * t / E),
    const   = function(t) rep(1, length(t))
  )
  out <- numeric(0)
  while (length(out) < n) {
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
run_method <- function(f) {
  res <- tryCatch(
    suppressWarnings({
      setTimeLimit(elapsed = time_limit, transient = TRUE)
      on.exit(setTimeLimit(elapsed = Inf), add = TRUE)
      f()
    }),
    error = function(e) NULL)
  if (is.null(res)) return(c(estimate = NA_real_, lower = NA_real_, upper = NA_real_))
  c(estimate = res$estimate, lower = res$lower, upper = res$upper)
}

# --- part k -------------------------------------------------------------------
k_methods <- list(
  ole                   = function(sd, k) ole(sd, alpha = alpha, k = k),
  hall_wang1999         = function(sd, k) hall_wang1999(sd, alpha = alpha, k = k),
  fraga_alves_neves2014 = function(sd, k) fraga_alves_neves2014(sd, k = k),
  dekkers1989           = function(sd, k) dekkers1989(sd, k = k, alpha = alpha)
)
k_task <- function(task) {
  set.seed(seed_base + task$id)
  sd <- make_record(task$n, task$scenario, E)
  n_events <- sum(sd$count); n_dist <- nrow(sd)
  k_sel <- tryCatch(suppressWarnings(
    ede_select_k(sd, nsim = 199L, alpha = 0.1)$k), error = function(e) NA_integer_)
  rules <- list(`5` = 5L, `10` = 10L, `20` = 20L, all = NA_integer_, sel = k_sel)
  rows <- list()
  for (m in names(k_methods)) {
    avail <- if (m == "ole") n_events else n_dist
    cache <- list()
    for (rule in names(rules)) {
      k <- rules[[rule]]
      if (rule == "all") k <- avail
      if (is.na(k)) {
        r <- c(estimate = NA_real_, lower = NA_real_, upper = NA_real_); ku <- NA_integer_
      } else {
        ku <- as.integer(min(k, avail)); key <- as.character(ku)
        if (is.null(cache[[key]]))
          cache[[key]] <- run_method(function() k_methods[[m]](sd, ku))
        r <- cache[[key]]
      }
      rows[[length(rows) + 1L]] <- data.frame(
        part = "k", scenario = task$scenario, n = task$n, rep = task$rep,
        method = m, k_rule = rule, k_used = ku,
        estimate = r[["estimate"]], lower = r[["lower"]], upper = r[["upper"]],
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

# --- part extant ----------------------------------------------------------------
test_year_e <- 101
extant_methods <- function(sd) list(
  ole               = function() ole(sd, alpha = alpha),
  robson1964        = function() robson1964(sd, alpha = alpha),
  strauss1989       = function() strauss1989(sd, alpha = alpha),
  solow1993         = function() solow1993(sd, alpha = alpha, test_year = test_year_e),
  solow1993b        = function() solow1993b(sd, alpha = alpha, test_year = test_year_e),
  solow2005         = function() solow2005(sd, alpha = alpha, test_year = test_year_e),
  mcinerny2006      = function() mcinerny2006(sd, alpha = alpha, test_year = test_year_e),
  burgman1995       = function() burgman1995(sd, alpha = alpha, test_year = test_year_e),
  solow_roberts2003 = function() solow_roberts2003(sd, alpha = alpha, test_year = test_year_e),
  jaric2010         = function() jaric2010(sd, alpha = alpha, test_year = test_year_e),
  hall_wang1999     = function() hall_wang1999(sd, alpha = alpha)
)
extant_task <- function(task) {
  set.seed(seed_base + task$id)
  sd <- make_record(task$n, task$scenario, E)
  calls <- extant_methods(sd)
  rows <- lapply(names(calls), function(m) {
    r <- run_method(calls[[m]])
    data.frame(part = "extant", scenario = task$scenario, n = task$n, rep = task$rep,
               method = m, k_rule = NA_character_, k_used = NA_integer_,
               estimate = r[["estimate"]], lower = r[["lower"]], upper = r[["upper"]],
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# --- cluster ----------------------------------------------------------------------
pkg_root <- normalizePath(".")
cl <- parallel::makeCluster(cores)
on.exit(parallel::stopCluster(cl), add = TRUE)
parallel::clusterExport(cl, c("E", "alpha", "time_limit", "seed_base", "test_year_e",
                              "draw_times", "make_record", "run_method", "k_methods",
                              "k_task", "extant_methods", "extant_task", "pkg_root"))
invisible(parallel::clusterEvalQ(
  cl, { suppressMessages(devtools::load_all(pkg_root, quiet = TRUE)); NULL }))

run_cells <- function(part_name, scenarios, sizes, reps, task_fun, id_offset) {
  for (sc in scenarios) for (n in sizes) {
    f <- file.path(outdir, sprintf("%s_%s_n%d.rds", part_name, sc, n))
    if (file.exists(f)) { cat("skip", basename(f), "\n"); next }
    t0 <- Sys.time()
    id0 <- id_offset + match(sc, scenarios) * 1e6 + n * 1e3
    tasks <- lapply(seq_len(reps), function(i)
      list(id = id0 + i, scenario = sc, n = n, rep = i))
    res <- do.call(rbind, parallel::parLapplyLB(cl, tasks, task_fun))
    saveRDS(res, f)
    cat(sprintf("%-7s %-8s n=%-3d %d reps  %.1f min\n", part_name, sc, n, reps,
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}

if (part %in% c("k", "both"))
  run_cells("k", c("uniform", "beta2", "beta4", "eff_dec", "eff_inc"),
            c(10L, 20L, 50L), reps_k, k_task, 1e8)
if (part %in% c("extant", "both"))
  run_cells("extant", c("const", "decl2", "decl3"),
            c(5L, 10L, 20L, 50L), reps_e, extant_task, 2e8)

# --- summaries ----------------------------------------------------------------------
read_part <- function(p) {
  f <- list.files(outdir, pattern = sprintf("^%s_.*_n[0-9]+\\.rds$", p), full.names = TRUE)
  if (length(f)) do.call(rbind, lapply(f, readRDS)) else NULL
}
kk <- read_part("k")
if (!is.null(kk)) {
  saveRDS(kk, file.path(outdir, "simulation_k_raw.rds"))
  cell <- function(d) {
    ok <- !is.na(d$estimate) | !is.na(d$upper)
    has_int <- !(d$method[1] %in% c("fraga_alves_neves2014", "dekkers1989"))
    cov <- ok & (is.na(d$lower) | d$lower <= E) & (is.na(d$upper) | d$upper >= E)
    fin <- ok & is.finite(d$upper)
    data.frame(scenario = d$scenario[1], n = d$n[1], method = d$method[1],
               k_rule = d$k_rule[1], reps = nrow(d), failure = mean(!ok),
               k_used_median = median(d$k_used, na.rm = TRUE),
               bias_median = if (any(!is.na(d$estimate))) median(d$estimate - E, na.rm = TRUE) else NA_real_,
               mae_median = if (any(!is.na(d$estimate))) median(abs(d$estimate - E), na.rm = TRUE) else NA_real_,
               coverage_all = if (has_int) sum(cov) / nrow(d) else NA_real_,
               share_finite = if (has_int && any(ok)) sum(fin) / sum(ok) else NA_real_,
               stringsAsFactors = FALSE)
  }
  s <- do.call(rbind, lapply(split(kk, list(kk$scenario, kk$n, kk$method, kk$k_rule), drop = TRUE), cell))
  rownames(s) <- NULL
  write.csv(s, file.path(outdir, "simulation_k_summary.csv"), row.names = FALSE)
  cat("k summary:", nrow(s), "rows\n")
}
ee <- read_part("extant")
if (!is.null(ee)) {
  saveRDS(ee, file.path(outdir, "simulation_extant_raw.rds"))
  cell <- function(d) {
    ok <- !is.na(d$estimate) | !is.na(d$upper)
    alarm <- !is.na(d$upper) & d$upper <= test_year_e
    data.frame(scenario = d$scenario[1], n = d$n[1], method = d$method[1],
               reps = nrow(d), failure = mean(!ok), false_alarm = mean(alarm),
               stringsAsFactors = FALSE)
  }
  s <- do.call(rbind, lapply(split(ee, list(ee$scenario, ee$n, ee$method), drop = TRUE), cell))
  rownames(s) <- NULL
  write.csv(s, file.path(outdir, "simulation_extant_summary.csv"), row.names = FALSE)
  cat("extant summary:", nrow(s), "rows\n")
}
cat("Done. Files in", normalizePath(outdir), "\n")
