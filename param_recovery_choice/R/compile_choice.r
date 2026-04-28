#!/usr/bin/env Rscript
setwd("~/reliable_info/param_recovery_choice")
renv::load()
options(renv.config.install.staged = FALSE)
renv::install(c("dplyr", "tibble", "ggplot2", "patchwork", "stan-dev/cmdstanr"))
renv::snapshot()
library(cmdstanr)

model <- cmdstan_model(
  stan_file       = './stan/log_seq_basic.stan',
  force_recompile = TRUE,
  cpp_options     = list(stan_threads = TRUE),
  stanc_options   = list("O1")
)
cat(sprintf("Compiled: %s\n", model$exe_file()))
