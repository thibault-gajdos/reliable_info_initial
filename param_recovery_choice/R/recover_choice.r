rm(list = ls(all = TRUE))
setwd("~/reliable_info/param_recovery_choice")
renv::load()

library(cmdstanr)
library(dplyr)
library(tibble)

args <- commandArgs(trailingOnly = TRUE)
k    <- as.numeric(args[1])

# ============================================================================
# 0. LOAD PRE-COMPILED MODEL
# ============================================================================

model <- cmdstan_model(
  stan_file     = './stan/log_seq_basic.stan',
  cpp_options   = list(stan_threads = TRUE),
  stanc_options = list("O1")
)

# ============================================================================
# 1. LOAD STIMULUS DATA FROM ACTUAL EXPERIMENT
# ============================================================================

load(file = './data/data_reliability.rdata')

data <- data %>%
  mutate(choice = case_when(
    (ResponseButtonOrder == 1 & Response == 0) ~ 2,
    (ResponseButtonOrder == 1 & Response == 1) ~ 1,
    (ResponseButtonOrder == 0 & Response == 0) ~ 1,
    (ResponseButtonOrder == 0 & Response == 1) ~ 2
  )) %>%
  mutate_at(vars(starts_with("color")), ~ ifelse(. == "blue", 1, 2))

subjs  <- unique(data$ParticipantPrivateID)
N      <- length(subjs)
T_max  <- max(data$TrialNumber)
I      <- sum(grepl("^color_", names(data)))

d      <- data %>% group_by(ParticipantPrivateID) %>% summarise(t_subjs = n())
Tsubj  <- d$t_subjs

## Stimulus arrays
color_arr <- array(-1L,  c(N, T_max, I))
proba_arr <- array(-1.0, c(N, T_max, I))

for (n in seq_len(N)) {
  t_n       <- Tsubj[n]
  data_subj <- data %>% filter(ParticipantPrivateID == subjs[n])
  for (j in seq_len(t_n)) {
    for (i in seq_len(I)) {
      color_arr[n, j, i] <- data_subj[[paste0("color_", i)]][j]
      proba_arr[n, j, i] <- data_subj[[paste0("proba_", i)]][j] / 100
    }
  }
}

# ============================================================================
# 2. PARAMETER GRID  (from paper, 64 combinations)
# ============================================================================
#
# Stan transforms:
#   alpha = 20 * Phi(mu_pr[1] + sigma[1] * raw)     → alpha in [0, 20]
#   beta  = mu_pr[2] + sigma[2] * raw               → unbounded
#   w[i]  = 2  * Phi(mu_w[i]  + sigma_w[i] * raw)   → w in [0, 2]

alpha_range <- c(0.5, 1, 1.5, 2.5)
beta_range  <- c(-0.1, 0, 0.3, 0.5)

w_profiles <- matrix(
  c(1.0, 1.0, 1.0, 1.0, 1.0,
    0.5, 0.6, 0.7, 0.8, 0.9,
    1.5, 1.4, 1.3, 1.2, 1.0,
    1.0, 0.9, 0.8, 0.8, 0.9),
  nrow = 4, byrow = TRUE
)

parameters <- expand.grid(
  alpha_idx = seq_along(alpha_range),
  beta_idx  = seq_along(beta_range),
  w_profile = seq_len(nrow(w_profiles))
)

stopifnot(k >= 1, k <= nrow(parameters))

mu_alpha <- alpha_range[parameters$alpha_idx[k]]
mu_beta  <- beta_range[parameters$beta_idx[k]]
mu_w     <- w_profiles[parameters$w_profile[k], ]

## Convert to raw space to match Stan parameterization
mu_pr <- c(
  qnorm(mu_alpha / 20),    # alpha = 20 * Phi(mu_pr[1])
  mu_beta                   # beta  = mu_pr[2] (identity)
)
mu_w_pr <- qnorm(mu_w / 2) # w = 2 * Phi(mu_w_pr)

## Group-level SDs in raw space
sigma_pr <- c(0.4, 0.3)    # alpha, beta
sigma_w  <- rep(0.3, 5)

## True group-level values (interpretable) for recovery comparison
group_sim   <- c(mu_alpha, mu_beta, mu_w)
param_names <- c("alpha", "beta", "w1", "w2", "w3", "w4", "w5")

cat(sprintf("[k=%d] alpha=%.2f beta=%.2f w=[%s]\n",
            k, mu_alpha, mu_beta,
            paste(sprintf("%.2f", mu_w), collapse=",")))

# ============================================================================
# 3. SIMULATE RESPONSES (matching Stan's generative model exactly)
# ============================================================================

safe_logit <- function(p) {
  p <- pmax(pmin(p, 1 - 1e-6), 1e-6)
  log(p / (1 - p))
}

set.seed(42 + k)

## Draw individual parameters via Stan's non-centered parameterization
params_indiv <- matrix(NA_real_, nrow = N, ncol = 7,
                       dimnames = list(NULL, c("alpha","beta",
                                               "w1","w2","w3","w4","w5")))
for (n in seq_len(N)) {
  raw <- rnorm(7, 0, 1)
  params_indiv[n, 1] <- 20 * pnorm(mu_pr[1] + sigma_pr[1] * raw[1])   # alpha
  params_indiv[n, 2] <- mu_pr[2] + sigma_pr[2] * raw[2]                # beta
  for (j in 1:5)
    params_indiv[n, 2 + j] <- 2 * pnorm(mu_w_pr[j] + sigma_w[j] * raw[2 + j])  # w1-w5
}

## Simulate choices
choice_sim <- array(-1L, c(N, T_max))

for (n in seq_len(N)) {
  alpha_n <- params_indiv[n, 1]
  beta_n  <- params_indiv[n, 2]
  w_n     <- c(params_indiv[n, 3:7], 1.0)   # w6 = 1 fixed

  for (j in seq_len(Tsubj[n])) {
    evidence <- c(0.0, 0.0)
    for (s in seq_len(I)) {
      c_s <- color_arr[n, j, s]
      p_s <- proba_arr[n, j, s]
      if (c_s < 1 || c_s > 2 || p_s <= 0 || p_s >= 1) next
      m_s <- alpha_n * safe_logit(p_s) + beta_n
      evidence[c_s] <- evidence[c_s] + w_n[s] * m_s
    }
    ## softmax(evidence) — no temperature, matching Stan
    val <- exp(evidence - max(evidence))
    val <- val / sum(val)
    choice_sim[n, j] <- sample(1:2, 1, prob = val)
  }
}

# ============================================================================
# 4. FIT
# ============================================================================

dir.create('./results/choice', recursive = TRUE, showWarnings = FALSE)

out_file <- sprintf('./results/choice/recover_%d.rds', k)
if (file.exists(out_file)) {
  cat(sprintf("[k=%d] Output already exists — skipping.\n", k))
  quit(save = "no", status = 0)
}

fit <- model$sample(
  data = list(
    N      = N,
    T_max  = T_max,
    I      = I,
    Tsubj  = Tsubj,
    color  = color_arr,
    proba  = proba_arr,
    choice = choice_sim
  ),
  iter_sampling     = 3000,
  iter_warmup       = 2000,
  chains            = 4,
  parallel_chains   = 4,
  threads_per_chain = max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")) %/% 4),
  seed              = 12345 + k,
  adapt_delta       = 0.95,
  max_treedepth     = 12,
  refresh           = 500
)

# ============================================================================
# 5. SAVE
# ============================================================================

## Group-level recovery
group_param_names <- c('mu_alpha','mu_beta','mu_w1','mu_w2','mu_w3','mu_w4','mu_w5')
group_fitted <- fit$summary(variables = group_param_names)

## Individual-level recovery: alpha[n], beta[n], w[n,1..5]
indiv_alpha <- paste0("alpha[", seq_len(N), "]")
indiv_beta  <- paste0("beta[", seq_len(N), "]")
indiv_w     <- paste0("w[", rep(seq_len(N), each = 5), ",",
                       rep(1:5, times = N), "]")
indiv_vars   <- c(indiv_alpha, indiv_beta, indiv_w)
indiv_fitted <- fit$summary(variables = indiv_vars)

## Simulated individual values (matching order above)
indiv_sim_vec <- c(
  params_indiv[, 1],                # alpha for all N
  params_indiv[, 2],                # beta for all N
  as.vector(t(params_indiv[, 3:7])) # w[1,1]..w[1,5], w[2,1]..w[2,5], ...
)

## Diagnostics
diag_summary <- fit$diagnostic_summary(quiet = TRUE)

results <- list(
  k             = k,
  group_sim     = group_sim,
  group_fitted  = group_fitted,
  indiv_sim     = indiv_sim_vec,
  indiv_fitted  = indiv_fitted,
  param_names   = param_names,
  n_divergent   = sum(diag_summary$num_divergent),
  n_max_td      = sum(diag_summary$num_max_treedepth),
  max_rhat      = max(fit$summary()$rhat, na.rm = TRUE),
  min_ess       = min(fit$summary()$ess_bulk, na.rm = TRUE)
)

saveRDS(results, out_file)
cat(sprintf("[k=%d] Done. divergences=%d max_rhat=%.3f\n",
            k, results$n_divergent, results$max_rhat))
