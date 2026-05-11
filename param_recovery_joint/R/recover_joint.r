rm(list = ls(all = TRUE))
setwd("~/reliable_info/param_recovery_joint")
renv::load()

library(cmdstanr)
library(dplyr)
library(stringr)
library(tibble)

args <- commandArgs(trailingOnly = TRUE)
k    <- as.numeric(args[1])

# ============================================================================
# 0. LOAD PRE-COMPILED MODEL
# ============================================================================

model <- cmdstan_model(
  stan_file   = './stan/linear_choice_influence_proba.stan',
  cpp_options = list(stan_threads = TRUE),
  stanc_options = list("O1")
)

# ============================================================================
# 1. LOAD STIMULUS DATA FROM ACTUAL EXPERIMENT
# ============================================================================

load(file = './data/data_reliability_exp6.rdata')

## Recode choice, construct influence_sample indicator
data <- data %>%
  mutate(choice = case_when(
    (ResponseButtonOrder == 1 & Response == 0) ~ 2,
    (ResponseButtonOrder == 1 & Response == 1) ~ 1,
    (ResponseButtonOrder == 0 & Response == 0) ~ 1,
    (ResponseButtonOrder == 0 & Response == 1) ~ 2
  )) %>%
  mutate_at(vars(starts_with("color")), ~ ifelse(. == "blue", 1, 2)) %>%
  rowwise() %>%
  mutate(sample_number = sum(!is.na(c_across(starts_with("proba_"))))) %>%
  ungroup() %>%
  rename(influence_proba = SliderReliability, influence = SliderResponse) %>%
  mutate(influence = influence / 100 - 0.5) %>%
  mutate(across(proba_1:proba_6, ~ as.integer(.x == influence_proba), .names = "{.col}_inf")) %>%
  rename_with(~ str_replace(.x, "proba_(\\d+)_inf", "influence_\\1"))

subjs  <- unique(data$ParticipantPrivateID)
N      <- length(subjs)
T_max  <- max(data$TrialNumber)
I_max  <- max(data$sample_number)

d      <- data %>% group_by(ParticipantPrivateID) %>% summarise(t_subjs = n())
Tsubj  <- d$t_subjs

## Stimulus arrays
sample_arr      <- array(-1L,   c(N, T_max))
color_arr       <- array(-1L,   c(N, T_max, I_max))
proba_arr       <- array(-1.0,  c(N, T_max, I_max))
infl_sample_arr <- array(-1L,   c(N, T_max, I_max))

for (n in seq_len(N)) {
  t_n       <- Tsubj[n]
  data_subj <- data %>% filter(ParticipantPrivateID == subjs[n])

  for (j in seq_len(t_n)) {
    sample_arr[n, j] <- data_subj$sample_number[j]
    for (i in seq_len(data_subj$sample_number[j])) {
      color_arr[n, j, i]       <- data_subj[[paste0("color_", i)]][j]
      proba_arr[n, j, i]       <- data_subj[[paste0("proba_", i)]][j] / 100
      infl_sample_arr[n, j, i] <- data_subj[[paste0("influence_", i)]][j]
    }
  }
}

# ============================================================================
# 2. PARAMETER GRID
# ============================================================================
#
# Empirical estimates (for reference):
#   mu_w     ≈ (0.89, 0.90, 0.94, 0.95, 1.03)
#   mu_alpha ≈ 2.56,  mu_beta ≈ 0.28
#   mu_a_infl ≈ 0.17, mu_b_infl ≈ 0.13, mu_sigma_infl ≈ 0.14
#
# Grid: 4 alpha × 4 beta × 3 a_infl × 2 b_infl × 2 sigma_infl × 4 w = 768

alpha_range      <- c(1, 2, 3.5, 5)
beta_range       <- c(-0.2, 0, 0.3, 0.6)
a_infl_range     <- c(0, 0.15, 0.4)
b_infl_range     <- c(0, 0.15)
sigma_infl_range <- c(0.08, 0.15)

w_profiles <- matrix(
  c(1.0, 1.0,  1.0,  1.0,  1.0,    # uniform
    0.8, 0.85, 0.9,  0.95, 1.05,   # increasing (close to empirical)
    1.2, 1.1,  1.0,  0.9,  0.8,    # decreasing
    1.0, 0.8,  0.8,  0.8,  1.0),   # U-shape
  nrow = 4, byrow = TRUE
)

parameters <- expand.grid(
  alpha_idx      = seq_along(alpha_range),
  beta_idx       = seq_along(beta_range),
  a_infl_idx     = seq_along(a_infl_range),
  b_infl_idx     = seq_along(b_infl_range),
  sigma_infl_idx = seq_along(sigma_infl_range),
  w_profile      = seq_len(nrow(w_profiles))
)

stopifnot(k >= 1, k <= nrow(parameters))

mu_alpha      <- alpha_range[parameters$alpha_idx[k]]
mu_beta       <- beta_range[parameters$beta_idx[k]]
mu_a_infl     <- a_infl_range[parameters$a_infl_idx[k]]
mu_b_infl     <- b_infl_range[parameters$b_infl_idx[k]]
mu_sigma_infl <- sigma_infl_range[parameters$sigma_infl_idx[k]]
mu_w          <- w_profiles[parameters$w_profile[k], ]

## Group-level SDs in raw space (from empirical fit)
SIGMA_PR <- c(
  0.202,   # w1
  0.093,   # w2
  0.152,   # w3
  0.159,   # w4
  0.096,   # w5
  0.493,   # alpha
  0.413,   # beta
  0.308,   # a_infl
  0.078,   # b_infl
  0.335    # sigma_infl (on log scale)
)

## Convert to raw (probit/log) space to match Stan parameterization
mu_pr <- c(
  qnorm(mu_w / 2),           # w1-w5:       w = 2 * Phi(mu_pr)
  qnorm(mu_alpha / 10),      # alpha:   alpha = 10 * Phi(mu_pr)
  mu_beta,                    # beta:    identity
  mu_a_infl,                  # a_infl:  identity
  mu_b_infl,                  # b_infl:  identity
  log(mu_sigma_infl)          # sigma_infl: exp transform
)

## True group-level values (interpretable) for recovery comparison
group_sim   <- c(mu_w, mu_alpha, mu_beta, mu_a_infl, mu_b_infl, mu_sigma_infl)
param_names <- c("w1","w2","w3","w4","w5",
                 "alpha","beta","a_infl","b_infl","sigma_infl")

cat(sprintf("[k=%d] alpha=%.1f beta=%.1f a_infl=%.2f b_infl=%.2f sigma_infl=%.2f w=[%s]\n",
            k, mu_alpha, mu_beta, mu_a_infl, mu_b_infl, mu_sigma_infl,
            paste(sprintf("%.2f", mu_w), collapse=",")))

# ============================================================================
# 3. SIMULATE RESPONSES (matching Stan's generative model exactly)
# ============================================================================

safe_logit <- function(p) {
  p <- pmax(pmin(p, 1 - 1e-6), 1e-6)
  log(p / (1 - p))
}

get_weight <- function(s, w) if (s <= 5) w[s] else 1.0

simulate_trial <- function(j, n, params_n) {
  w            <- params_n[1:5]
  alpha_n      <- params_n[6]
  beta_n       <- params_n[7]
  a_infl_n     <- params_n[8]
  b_infl_n     <- params_n[9]
  sigma_infl_n <- params_n[10]

  S <- sample_arr[n, j]
  if (S < 1) return(list(choice = -1L, influence = -99.0))

  ev     <- c(0.0, 0.0)
  sub_ev <- 0.0

  for (s in seq_len(S)) {
    p <- proba_arr[n, j, s]
    c <- color_arr[n, j, s]
    if (c < 1 || c > 2 || p <= 0 || p >= 1) next

    m  <- alpha_n * safe_logit(p) + beta_n
    ws <- get_weight(s, w)
    ev[c] <- ev[c] + ws * m

    if (infl_sample_arr[n, j, s] == 1L) {
      d_col  <- if (c == 1L) 1.0 else -1.0
      sub_ev <- sub_ev + ws * m * d_col
    }
  }

  ev_diff <- ev[1] - ev[2]
  ch      <- if (runif(1) < plogis(ev_diff)) 1L else 2L

  ## Influence report (untruncated, matching Stan likelihood)
  total_ev_ch   <- if (ch == 1L) ev_diff else -ev_diff
  sub_ev_ch     <- if (ch == 1L) sub_ev  else -sub_ev
  non_sub_ev_ch <- total_ev_ch - sub_ev_ch
  delta_p       <- plogis(total_ev_ch) - plogis(non_sub_ev_ch)
  mu_infl       <- a_infl_n * delta_p + b_infl_n
  infl_rep      <- rnorm(1, mu_infl, sigma_infl_n)

  list(choice = ch, influence = infl_rep)
}

set.seed(42 + k)

## Draw individual parameters via Stan's non-centered parameterization
params_indiv <- matrix(NA_real_, nrow = N, ncol = 10,
                       dimnames = list(NULL, param_names))

for (n in seq_len(N)) {
  raw <- rnorm(10, 0, 1)
  for (j in 1:5)
    params_indiv[n, j] <- 2 * pnorm(mu_pr[j] + SIGMA_PR[j] * raw[j])
  params_indiv[n, 6]  <- 10 * pnorm(mu_pr[6]  + SIGMA_PR[6]  * raw[6])
  params_indiv[n, 7]  <- mu_pr[7]  + SIGMA_PR[7]  * raw[7]
  params_indiv[n, 8]  <- mu_pr[8]  + SIGMA_PR[8]  * raw[8]
  params_indiv[n, 9]  <- mu_pr[9]  + SIGMA_PR[9]  * raw[9]
  params_indiv[n, 10] <- exp(mu_pr[10] + SIGMA_PR[10] * raw[10])
}

## Simulate choices and influence reports
choice_sim    <- array(-1L,    c(N, T_max))
influence_sim <- array(-99.0,  c(N, T_max))

for (n in seq_len(N)) {
  for (j in seq_len(Tsubj[n])) {
    res <- simulate_trial(j, n, params_indiv[n, ])
    choice_sim[n, j]    <- res$choice
    influence_sim[n, j] <- res$influence
  }
}

# ============================================================================
# 4. FIT
# ============================================================================

dir.create('./results/joint', recursive = TRUE, showWarnings = FALSE)

out_file <- sprintf('./results/joint/recover_%d.rds', k)
if (file.exists(out_file)) {
  cat(sprintf("[k=%d] Output already exists — skipping.\n", k))
  quit(save = "no", status = 0)
}

fit <- model$sample(
  data = list(
    N                = N,
    T_max            = T_max,
    I_max            = I_max,
    Tsubj            = Tsubj,
    sample           = sample_arr,
    color            = color_arr,
    proba            = proba_arr,
    choice           = choice_sim,
    influence_sample = infl_sample_arr,
    influence        = influence_sim,
    grainsize        = 5
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
group_param_names <- c('mu_w1','mu_w2','mu_w3','mu_w4','mu_w5',
                       'mu_alpha','mu_beta','mu_a_infl','mu_b_infl','mu_sigma_infl')
group_fitted <- fit$summary(variables = group_param_names)

## Individual-level recovery
indiv_vars <- paste0("params[", rep(seq_len(N), each = 10), ",",
                     rep(1:10, times = N), "]")
indiv_fitted <- fit$summary(variables = indiv_vars)

## Diagnostics
diag_summary <- fit$diagnostic_summary(quiet = TRUE)

results <- list(
  k             = k,
  group_sim     = group_sim,
  group_fitted  = group_fitted,
  indiv_sim     = params_indiv,
  indiv_fitted  = indiv_fitted,
  param_names   = param_names,
  n_divergent   = sum(diag_summary$num_divergent),
  n_max_td      = sum(diag_summary$num_max_treedepth),
  max_rhat      = max(c(group_fitted$rhat, indiv_fitted$rhat), na.rm = TRUE),
  min_ess       = min(c(group_fitted$ess_bulk, indiv_fitted$ess_bulk), na.rm = TRUE)
)

saveRDS(results, sprintf('./results/joint/recover_%d.rds', k))
cat(sprintf("[k=%d] Done. divergences=%d max_rhat=%.3f\n",
            k, results$n_divergent, results$max_rhat))
