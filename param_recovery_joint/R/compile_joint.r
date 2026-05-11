#!/usr/bin/env Rscript
# ============================================================================
# Compile Stan model for parameter recovery
# ============================================================================
# Run ONCE before launching the SLURM array job:
#   Rscript compile_joint.r
#
# cmdstanr caches the binary next to the .stan file. Subsequent calls to
# cmdstan_model() with the same stan_file + options will reuse it.
# ============================================================================
setwd("~/reliable_info/param_recovery_joint")
##renv::restore()
renv::load()
library(cmdstanr)


stan_file <- "stan/linear_choice_influence_proba.stan"

if (!file.exists(stan_file))
  stop(sprintf("Stan file not found: %s", stan_file))

cat(sprintf("Compiling %s ...\n", stan_file))
model <- cmdstan_model(
  stan_file         = stan_file,
  force_recompile   = TRUE,
  cpp_options       = list(stan_threads = TRUE),
  stanc_options     = list("O1")
)

cat(sprintf("Compiled executable: %s\n", model$exe_file()))
cat("Done.\n")
