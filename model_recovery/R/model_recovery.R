#!/usr/bin/env Rscript

# ============================================================================
# Model Recovery: Simulate & Fit (cmdstanr)
# ============================================================================
# Simulates data from one generating model, fits all 4 model variants,
# saves fits and LOO objects. No model comparison — that's handled
# separately by compare_recovery.R.
#
# Job mapping:
#   --job_id indexes into the master grid of all valid combinations × reps.
#   The script maps this to (generating_model, param_combo, rep).
#
# Jitter strategy:
#   Grid values define the center. For each rep, group-level means are
#   jittered by small Gaussian noise (only for active dimensions).
#
# Hierarchical structure:
#   Individual parameters are generated through the same non-centered
#   parameterization as the Stan model (probit transforms + raw ~ N(0,1)).
#
# Seed strategy:
#   sim_seed  -> R RNG for jitter + data simulation (set.seed)
#   mcmc_seed -> passed explicitly to Stan $sample()
#   Both depend on (generating_model, param_combo, rep) for uniqueness.
# ============================================================================

setwd("~/reliable_info/model_recovery")
renv::load()

library(cmdstanr)
library(loo)
library(optparse)

# ============================================================================
# CLI ARGUMENTS
# ============================================================================

option_list <- list(
  make_option(c("-j", "--job_id"), type="integer", default=NULL,
              help="Master grid job ID (1 to 64*n_reps). Maps to (model, combo, rep)."),
  make_option(c("-g", "--generating_model"), type="integer", default=1,
              help="Generating model (1-4). Ignored if --job_id is set."),
  make_option(c("-c", "--param_combo"), type="integer", default=1,
              help="Parameter combination ID. Ignored if --job_id is set."),
  make_option(c("-r", "--rep"), type="integer", default=1,
              help="Repetition number (1 to n_reps). Ignored if --job_id is set."),
  make_option(c("--n_reps"), type="integer", default=10,
              help="Number of repetitions per combination (default: 10)"),
  make_option(c("-o", "--output_dir"), type="character", default="results",
              help="Base output directory (fits/ and loo/ subdirs created automatically)"),
  make_option(c("--n_iter"), type="integer", default=2000,
              help="MCMC iterations per chain (default: 2000, giving 1000 post-warmup)"),
  make_option(c("--n_chains"), type="integer", default=4,
              help="Number of MCMC chains"),
  make_option(c("--threads_per_chain"), type="integer", default=4,
              help="Threads per chain for within-chain parallelism (default: 4)"),
  make_option(c("--n_cores"), type="integer", default=8,
              help="Total cores available. Used to compute parallel_chains."),
  make_option(c("--adapt_delta"), type="double", default=0.95,
              help="Stan adapt_delta (default: 0.95)"),
  make_option(c("--max_treedepth"), type="integer", default=12,
              help="Stan max_treedepth (default: 12)"),
  make_option(c("--stan_file"), type="character", default="stan/linear_choice_flexible.stan",
              help="Stan source file"),
  make_option(c("--mcmc_seed"), type="integer", default=NULL,
              help="Override MCMC seed (default: derived from job identity)")
)

opt_parser <- OptionParser(option_list=option_list)
opt        <- parse_args(opt_parser)

N_REPS            <- opt$n_reps
threads_per_chain <- opt$threads_per_chain
parallel_chains   <- min(opt$n_chains, max(1L, opt$n_cores %/% threads_per_chain))

# ============================================================================
# CONSTANTS
# ============================================================================

N_SUBJECTS <- 24
N_TRIALS   <- 300
I_MAX      <- 6
GRAINSIZE  <- 5

JITTER_SD_A <- 0.15
JITTER_SD_B <- 0.05
JITTER_SD_W <- 0.05

fits_dir <- file.path(opt$output_dir, "fits")
loo_dir  <- file.path(opt$output_dir, "loo")
dir.create(fits_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(loo_dir,  showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# MODEL CONFIGURATION
# ============================================================================

model_configs <- data.frame(
  model_id       = 1:4,
  name           = c("Full", "No_Weights", "No_Distortion", "Null"),
  fix_weights    = c(0L, 1L, 0L, 1L),
  fix_distortion = c(0L, 0L, 1L, 1L),
  stringsAsFactors = FALSE
)

# ============================================================================
# PARAMETER GRID
# ============================================================================

a_values <- c(0.5, 1.0, 1.5, 2.5)
b_values <- c(-0.1, 0.0, 0.3, 0.5)
w_sets   <- list(
  c(1, 1, 1, 1, 1),
  c(0.5, 0.6, 0.7, 0.8, 0.9),
  c(1.5, 1.4, 1.3, 1.2, 1.0),
  c(1, 0.9, 0.8, 0.8, 0.9)
)

# ============================================================================
# HIERARCHICAL STRUCTURE (must match Stan model)
# ============================================================================
#   w[j]  = 2  * Phi(mu_pr[j] + sigma_pr[j] * raw[j])   bounded [0, 2]
#   alpha = 10 * Phi(mu_pr[6] + sigma_pr[6] * raw[6])    bounded [0, 10]
#   beta  = mu_pr[7] + sigma_pr[7] * raw[7]              unbounded

SIGMA_PR_W <- 0.5
SIGMA_PR_A <- 0.5
SIGMA_PR_B <- 0.5

inv_w     <- function(w) qnorm(pmin(pmax(w, 0.001), 1.999) / 2)
inv_alpha <- function(a) qnorm(pmin(pmax(a, 0.001), 9.999) / 10)
fwd_w     <- function(mu, sig, raw) 2  * pnorm(mu + sig * raw)
fwd_alpha <- function(mu, sig, raw) 10 * pnorm(mu + sig * raw)
fwd_beta  <- function(mu, sig, raw) mu + sig * raw

# ============================================================================
# DEGENERATE COMBINATION FILTERING
# ============================================================================

full_grid <- expand.grid(a_idx = 1:4, b_idx = 1:4, w_idx = 1:4)
full_grid$degen_weights    <- full_grid$w_idx == 1
full_grid$degen_distortion <- (full_grid$a_idx == 2) & (full_grid$b_idx == 2)

build_valid_grid <- function(gen_model) {
  switch(gen_model,
    `1` = full_grid[!full_grid$degen_weights & !full_grid$degen_distortion, ],
    `2` = {
      g <- unique(full_grid[, c("a_idx", "b_idx")])
      g$w_idx <- NA_integer_; g$degen_weights <- FALSE
      g$degen_distortion <- (g$a_idx == 2) & (g$b_idx == 2)
      g[!g$degen_distortion, ]
    },
    `3` = {
      g <- data.frame(w_idx = 1:4)
      g$a_idx <- NA_integer_; g$b_idx <- NA_integer_
      g$degen_weights <- g$w_idx == 1; g$degen_distortion <- FALSE
      g[!g$degen_weights, ]
    },
    `4` = data.frame(a_idx = NA_integer_, b_idx = NA_integer_, w_idx = NA_integer_,
                     degen_weights = FALSE, degen_distortion = FALSE)
  )
}

# ============================================================================
# MASTER GRID (combos × reps)
# ============================================================================

combo_grid <- do.call(rbind, lapply(1:4, function(m) {
  vg <- build_valid_grid(m); rownames(vg) <- NULL
  vg$generating_model <- m; vg$local_combo <- seq_len(nrow(vg)); vg
}))
rownames(combo_grid) <- NULL

master_grid <- do.call(rbind, lapply(seq_len(N_REPS), function(r) {
  g <- combo_grid; g$rep <- r; g
}))
rownames(master_grid) <- NULL
master_grid$job_id <- seq_len(nrow(master_grid))
n_total_jobs <- nrow(master_grid)

cat(sprintf("Master grid: %d combos x %d reps = %d total jobs\n",
            nrow(combo_grid), N_REPS, n_total_jobs))

# ============================================================================
# RESOLVE JOB IDENTITY
# ============================================================================

if (!is.null(opt$job_id)) {
  if (!opt$job_id %in% master_grid$job_id)
    stop(sprintf("--job_id must be 1-%d. Got %d.", n_total_jobs, opt$job_id))
  row <- master_grid[master_grid$job_id == opt$job_id, ]
  opt$generating_model <- row$generating_model
  opt$param_combo      <- row$local_combo
  opt$rep              <- row$rep
} else {
  if (!opt$generating_model %in% 1:4) stop("--generating_model must be 1-4.")
  vg <- build_valid_grid(opt$generating_model)
  if (!opt$param_combo %in% seq_len(nrow(vg)))
    stop(sprintf("--param_combo out of range for model %d.", opt$generating_model))
  if (!opt$rep %in% seq_len(N_REPS))
    stop(sprintf("--rep must be 1-%d.", N_REPS))
}

gen_cfg    <- model_configs[opt$generating_model, ]
model_name <- gen_cfg$name
valid_grid <- build_valid_grid(opt$generating_model)
n_valid    <- nrow(valid_grid)
combo      <- valid_grid[opt$param_combo, ]

# ============================================================================
# SKIP IF ALREADY COMPLETED
# ============================================================================
# A job is complete if all 4 LOO files exist

all_loo_exist <- all(sapply(1:4, function(m) {
  file.exists(file.path(loo_dir,
    sprintf("loo_gen%d_combo%03d_rep%02d_fitted%d.rds",
            opt$generating_model, opt$param_combo, opt$rep, m)))
}))

if (all_loo_exist) {
  cat(sprintf("All 4 LOO files already exist for gen%d_combo%03d_rep%02d. Skipping.\n",
              opt$generating_model, opt$param_combo, opt$rep))
  quit(status = 0)
}

# ============================================================================
# SEEDS
# ============================================================================

sim_seed  <- (opt$generating_model * 10000L + opt$param_combo * 100L + opt$rep) * 9973L
mcmc_seed <- if (!is.null(opt$mcmc_seed)) opt$mcmc_seed else
             (opt$generating_model * 10000L + opt$param_combo * 100L + opt$rep) * 6271L
set.seed(sim_seed)

# ============================================================================
# JITTER GROUP-LEVEL MEANS
# ============================================================================

a_center <- if (!is.na(combo$a_idx)) a_values[combo$a_idx] else 1.0
b_center <- if (!is.na(combo$b_idx)) b_values[combo$b_idx] else 0.0
w_center <- if (!is.na(combo$w_idx)) w_sets[[combo$w_idx]] else rep(1.0, 5)

a_mean  <- a_center
b_mean  <- b_center
w_means <- w_center

if (!is.na(combo$a_idx)) a_mean  <- a_center + rnorm(1, 0, JITTER_SD_A)
if (!is.na(combo$b_idx)) b_mean  <- b_center + rnorm(1, 0, JITTER_SD_B)
if (!is.na(combo$w_idx)) w_means <- w_center + rnorm(5, 0, JITTER_SD_W)

a_mean  <- max(0.05, min(9.95, a_mean))
w_means <- pmax(0.05, pmin(1.95, w_means))

cat("=== Model Recovery: Simulate & Fit ===\n")
cat(sprintf("Generating model : %d (%s)\n", opt$generating_model, model_name))
cat(sprintf("Combo: %d/%d | Rep: %d/%d | Job: %s/%d\n",
            opt$param_combo, n_valid, opt$rep, N_REPS,
            ifelse(is.null(opt$job_id), "manual", as.character(opt$job_id)), n_total_jobs))
cat(sprintf("MCMC: %d chains x %d iter | threads=%d | parallel=%d\n",
            opt$n_chains, opt$n_iter, threads_per_chain, parallel_chains))
cat(sprintf("Seeds: sim=%d, mcmc=%d\n", sim_seed, mcmc_seed))
cat(sprintf("a=%.3f (center=%.2f) | b=%.3f (center=%.2f)\n",
            a_mean, a_center, b_mean, b_center))
cat(sprintf("w=[%s]\n", paste(sprintf("%.3f", w_means), collapse=", ")))
cat("\n")

# ============================================================================
# GENERATE INDIVIDUAL PARAMETERS (matching Stan hierarchy)
# ============================================================================

mu_pr    <- numeric(7)
sigma_pr <- numeric(7)
for (j in 1:5) { mu_pr[j] <- inv_w(w_means[j]); sigma_pr[j] <- SIGMA_PR_W }
mu_pr[6] <- inv_alpha(a_mean); sigma_pr[6] <- SIGMA_PR_A
mu_pr[7] <- b_mean;            sigma_pr[7] <- SIGMA_PR_B

param_raw <- matrix(rnorm(N_SUBJECTS * 7), nrow = N_SUBJECTS, ncol = 7)

w_indiv <- matrix(NA_real_, nrow = N_SUBJECTS, ncol = 5)
for (j in 1:5) w_indiv[, j] <- fwd_w(mu_pr[j], sigma_pr[j], param_raw[, j])
a_indiv <- fwd_alpha(mu_pr[6], sigma_pr[6], param_raw[, 6])
b_indiv <- fwd_beta(mu_pr[7], sigma_pr[7], param_raw[, 7])

cat(sprintf("alpha: mean=%.3f sd=%.3f | beta: mean=%.3f sd=%.3f\n",
            mean(a_indiv), sd(a_indiv), mean(b_indiv), sd(b_indiv)))

# ============================================================================
# SIMULATE DATA
# ============================================================================

cat("Simulating data...\n")

color_arr  <- array(-1L, dim = c(N_SUBJECTS, N_TRIALS, I_MAX))
proba_arr  <- array(-1,  dim = c(N_SUBJECTS, N_TRIALS, I_MAX))
sample_cnt <- matrix(I_MAX, nrow = N_SUBJECTS, ncol = N_TRIALS)
choice_mat <- matrix(-1L,   nrow = N_SUBJECTS, ncol = N_TRIALS)

map_prob   <- function(p, a, b, fix) if (fix == 1L) qlogis(p) else a * qlogis(p) + b
get_weight <- function(s, w, fix) if (fix == 1L) 1.0 else if (s <= 5L) w[s] else 1.0

for (n in seq_len(N_SUBJECTS)) {
  for (t in seq_len(N_TRIALS)) {
    color_arr[n, t, ] <- base::sample(1:2, I_MAX, replace = TRUE)
    proba_arr[n, t, ] <- runif(I_MAX, 0.35, 0.65)
    evidence <- c(0.0, 0.0)
    for (s in seq_len(I_MAX)) {
      c_s <- color_arr[n, t, s]; p_s <- proba_arr[n, t, s]
      m_s <- map_prob(p_s, a_indiv[n], b_indiv[n], gen_cfg$fix_distortion)
      w_s <- get_weight(s, w_indiv[n, ], gen_cfg$fix_weights)
      evidence[c_s] <- evidence[c_s] + w_s * m_s
    }
    choice_mat[n, t] <- ifelse(runif(1) < plogis(evidence[1] - evidence[2]), 1L, 2L)
  }
}

cat(sprintf("Choice distribution: %.1f%% blue\n\n", 100 * mean(choice_mat == 1L)))

# ============================================================================
# STAN DATA
# ============================================================================

stan_data_base <- list(
  N = N_SUBJECTS, T_max = N_TRIALS, I_max = I_MAX,
  Tsubj = rep(N_TRIALS, N_SUBJECTS), sample = sample_cnt,
  color = color_arr, proba = proba_arr, choice = choice_mat,
  grainsize = GRAINSIZE
)

# ============================================================================
# LOAD COMPILED MODEL
# ============================================================================

cat(sprintf("Loading Stan model: %s\n", opt$stan_file))
compiled_model <- cmdstan_model(
  stan_file = opt$stan_file,
  cpp_options = list(stan_threads = TRUE),
  stanc_options = list("O1"),
  compile_model_methods = TRUE
)
cat("\n")

# ============================================================================
# FIT ALL 4 MODEL VARIANTS
# ============================================================================

for (i in seq_len(nrow(model_configs))) {
  cfg <- model_configs[i, ]

  # Skip if this specific LOO file already exists
  loo_file <- file.path(loo_dir,
    sprintf("loo_gen%d_combo%03d_rep%02d_fitted%d.rds",
            opt$generating_model, opt$param_combo, opt$rep, cfg$model_id))
  if (file.exists(loo_file)) {
    cat(sprintf("--- Model %d (%s): LOO exists, skipping ---\n\n", cfg$model_id, cfg$name))
    next
  }

  cat(sprintf("--- Fitting Model %d: %s ---\n", cfg$model_id, cfg$name))

  stan_data_i <- c(stan_data_base, list(
    fix_sequential_weights     = cfg$fix_weights,
    fix_probability_distortion = cfg$fix_distortion
  ))

  fit <- tryCatch(
    compiled_model$sample(
      data = stan_data_i, chains = opt$n_chains,
      iter_sampling = opt$n_iter %/% 2L, iter_warmup = opt$n_iter %/% 2L,
      parallel_chains = parallel_chains, threads_per_chain = threads_per_chain,
      seed = mcmc_seed, refresh = 100, show_messages = FALSE,
      adapt_delta = opt$adapt_delta, max_treedepth = opt$max_treedepth
    ),
    error = function(e) { cat(sprintf("  ERROR: %s\n", e$message)); NULL }
  )

  if (is.null(fit)) { cat("  Skipping (fit failed).\n\n"); next }

  # Diagnostics
  diag <- tryCatch(fit$diagnostic_summary(quiet = TRUE), error = function(e) NULL)
  if (!is.null(diag)) {
    if (sum(diag$num_divergent) > 0)
      cat(sprintf("  WARNING: %d divergent transitions\n", sum(diag$num_divergent)))
    if (sum(diag$num_max_treedepth) > 0)
      cat(sprintf("  WARNING: %d max treedepth hits\n", sum(diag$num_max_treedepth)))
  }

  # Save fit
  fit_file <- file.path(fits_dir,
    sprintf("fit_gen%d_combo%03d_rep%02d_fitted%d.rds",
            opt$generating_model, opt$param_combo, opt$rep, cfg$model_id))
  fit$save_object(file = fit_file)
  cat(sprintf("  Fit saved: %s\n", fit_file))

  # LOO
  loo_obj <- tryCatch({
    log_lik_draws <- fit$draws("log_lik", format = "draws_array")
    r_eff <- relative_eff(exp(log_lik_draws))
    loo(log_lik_draws, r_eff = r_eff)
  }, error = function(e) {
    cat(sprintf("  LOO failed: %s\n", e$message)); NULL
  })

  if (!is.null(loo_obj)) {
    saveRDS(loo_obj, loo_file)
    cat(sprintf("  LOO saved | ELPD=%.2f (SE=%.2f) | p_loo=%.2f\n",
                loo_obj$estimates["elpd_loo", "Estimate"],
                loo_obj$estimates["elpd_loo", "SE"],
                loo_obj$estimates["p_loo", "Estimate"]))
    k_vals <- loo_obj$diagnostics$pareto_k
    cat(sprintf("  Pareto-k: %d>0.7, %d in (0.5,0.7]\n",
                sum(k_vals > 0.7, na.rm=TRUE), sum(k_vals > 0.5 & k_vals <= 0.7, na.rm=TRUE)))
  }

  cat("\n")
  rm(fit); gc()
}

cat("Done!\n")
