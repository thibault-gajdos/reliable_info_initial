#!/usr/bin/env Rscript

# ============================================================================
# Model Recovery: Compare & Analyze
# ============================================================================
# Reads all LOO files produced by model_recovery_v2.R, performs pairwise
# model comparison for each (generating_model, combo, rep), builds the
# confusion matrix, and produces plots + LaTeX table.
#
# Run after all fitting jobs have completed:
#   Rscript R/compare_recovery.R --loo_dir results/loo --output_dir results
# ============================================================================

setwd("~/reliable_info/model_recovery")
renv::load()

library(loo)
library(ggplot2)
library(dplyr)
library(tidyr)
library(optparse)

# ============================================================================
# CLI ARGUMENTS
# ============================================================================

option_list <- list(
  make_option("--loo_dir", type="character", default="results/loo",
              help="Directory containing LOO .rds files"),
  make_option("--output_dir", type="character", default="results",
              help="Directory for output files"),
  make_option("--n_reps", type="integer", default=10,
              help="Number of reps per combo (must match fitting script)")
)

opt <- parse_args(OptionParser(option_list=option_list))

# ============================================================================
# MODEL CONFIGURATION (must match model_recovery_v2.R)
# ============================================================================

model_configs <- data.frame(
  model_id       = 1:4,
  name           = c("Full", "No_Weights", "No_Distortion", "Null"),
  fix_weights    = c(0L, 1L, 0L, 1L),
  fix_distortion = c(0L, 0L, 1L, 1L),
  stringsAsFactors = FALSE
)
model_names <- model_configs$name
n_models    <- nrow(model_configs)

# ============================================================================
# DISCOVER ALL COMPLETED JOBS
# ============================================================================

cat("=== Model Recovery: Compare & Analyze ===\n")
cat(sprintf("LOO directory: %s\n\n", opt$loo_dir))

loo_files <- list.files(opt$loo_dir, pattern = "^loo_gen\\d+_combo\\d+_rep\\d+_fitted\\d+\\.rds$",
                        full.names = TRUE)
cat(sprintf("Found %d LOO files\n", length(loo_files)))

# Parse filenames to identify jobs
parse_loo_filename <- function(f) {
  m <- regmatches(basename(f), regexec("loo_gen(\\d+)_combo(\\d+)_rep(\\d+)_fitted(\\d+)\\.rds", basename(f)))[[1]]
  if (length(m) < 5) return(NULL)
  data.frame(gen = as.integer(m[2]), combo = as.integer(m[3]),
             rep = as.integer(m[4]), fitted = as.integer(m[5]),
             file = f, stringsAsFactors = FALSE)
}

file_info <- do.call(rbind, lapply(loo_files, parse_loo_filename))
if (is.null(file_info) || nrow(file_info) == 0)
  stop("No valid LOO files found!")

# Find jobs with all 4 fitted models
job_keys <- file_info %>%
  group_by(gen, combo, rep) %>%
  summarise(n_fitted = n_distinct(fitted), .groups = "drop")

complete_jobs <- job_keys %>% filter(n_fitted == 4)
incomplete_jobs <- job_keys %>% filter(n_fitted < 4)

cat(sprintf("Complete jobs (4/4 models): %d\n", nrow(complete_jobs)))
if (nrow(incomplete_jobs) > 0)
  cat(sprintf("Incomplete jobs (<4 models): %d (skipped)\n", nrow(incomplete_jobs)))
cat("\n")

if (nrow(complete_jobs) == 0) stop("No complete jobs found!")

# ============================================================================
# COMPARE MODELS FOR EACH JOB
# ============================================================================

cat("Running loo_compare for each job...\n")

comparison_results <- list()

for (r in seq_len(nrow(complete_jobs))) {
  job <- complete_jobs[r, ]

  # Load all 4 LOO objects for this job
  loo_list <- setNames(
    lapply(1:4, function(m) {
      f <- file_info %>%
        filter(gen == job$gen, combo == job$combo, rep == job$rep, fitted == m)
      if (nrow(f) == 0) return(NULL)
      tryCatch(readRDS(f$file[1]), error = function(e) NULL)
    }),
    model_names
  )

  # Skip if any failed to load
  loaded <- !sapply(loo_list, is.null)
  if (sum(loaded) < 2) next

  # Compare
  comp <- tryCatch(
    loo_compare(loo_list[loaded]),
    error = function(e) {
      cat(sprintf("  WARNING: loo_compare failed for gen%d_combo%03d_rep%02d: %s\n",
                  job$gen, job$combo, job$rep, e$message))
      NULL
    }
  )

  if (is.null(comp)) next

  winning_name <- rownames(comp)[1]
  winning_id   <- model_configs$model_id[model_configs$name == winning_name]

  # ELPD values for all fitted models
  elpd_vals <- sapply(loo_list[loaded], function(x) x$estimates["elpd_loo", "Estimate"])

  comparison_results[[length(comparison_results) + 1]] <- list(
    gen         = job$gen,
    combo       = job$combo,
    rep         = job$rep,
    winning_id  = winning_id,
    winning_name = winning_name,
    correct     = (winning_id == job$gen),
    elpd_values = elpd_vals,
    comparison  = comp
  )

  if (r %% 50 == 0) cat(sprintf("  Processed %d/%d jobs\n", r, nrow(complete_jobs)))
}

cat(sprintf("Successfully compared: %d jobs\n\n", length(comparison_results)))

if (length(comparison_results) == 0) stop("No successful comparisons!")

# ============================================================================
# BUILD DATA FRAME
# ============================================================================

recovery_df <- do.call(rbind, lapply(comparison_results, function(r) {
  data.frame(
    generating_model = r$gen,
    gen_model_name   = model_names[r$gen],
    param_combo      = r$combo,
    rep              = r$rep,
    winning_model    = r$winning_id,
    winning_name     = r$winning_name,
    correct          = r$correct,
    stringsAsFactors = FALSE
  )
}))

# ============================================================================
# CONFUSION MATRIX
# ============================================================================

confusion_matrix <- matrix(0L, nrow = n_models, ncol = n_models,
                            dimnames = list(model_names, model_names))

for (i in seq_len(nrow(recovery_df))) {
  g <- recovery_df$generating_model[i]
  w <- recovery_df$winning_model[i]
  if (g %in% 1:4 && w %in% 1:4)
    confusion_matrix[g, w] <- confusion_matrix[g, w] + 1L
}

row_totals       <- rowSums(confusion_matrix)
confusion_prop   <- prop.table(confusion_matrix, margin = 1)
inversion_matrix <- prop.table(confusion_matrix, margin = 2)
overall_accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)

cat("Confusion Matrix (raw counts):\n")
cat("Rows = Generating model, Columns = Selected model\n\n")
print(confusion_matrix)
cat("\n")

cat("Confusion Matrix (proportions):\n")
print(round(confusion_prop, 3))
cat("\n")

cat(sprintf("Overall recovery accuracy: %.1f%%\n", 100 * overall_accuracy))
for (i in seq_len(n_models)) {
  cat(sprintf("  %s: %.1f%% (%d/%d)\n",
              model_names[i],
              100 * confusion_matrix[i, i] / max(row_totals[i], 1),
              confusion_matrix[i, i], row_totals[i]))
}
cat("\n")

cat("Inversion Matrix: P(Generated by X | Selected Y)\n")
print(round(inversion_matrix, 3))
cat("\n")

# ============================================================================
# ACCURACY BY MODEL (with CI)
# ============================================================================

accuracy_by_model <- recovery_df %>%
  group_by(generating_model, gen_model_name) %>%
  summarise(
    accuracy  = mean(correct),
    n_correct = sum(correct),
    n_total   = n(),
    ci_lower  = binom.test(sum(correct), n())$conf.int[1],
    ci_upper  = binom.test(sum(correct), n())$conf.int[2],
    .groups   = "drop"
  )

cat("Recovery accuracy by model (95% CI):\n")
for (i in seq_len(nrow(accuracy_by_model))) {
  r <- accuracy_by_model[i, ]
  cat(sprintf("  %-15s %.1f%% [%.1f%%, %.1f%%]  (%d/%d)\n",
              r$gen_model_name,
              100 * r$accuracy, 100 * r$ci_lower, 100 * r$ci_upper,
              r$n_correct, r$n_total))
}
cat("\n")

# ============================================================================
# SAVE SUMMARY
# ============================================================================

summary_file <- file.path(opt$output_dir, "recovery_summary.rds")
saveRDS(list(
  confusion_matrix      = confusion_matrix,
  confusion_proportions = confusion_prop,
  inversion_matrix      = inversion_matrix,
  overall_accuracy      = overall_accuracy,
  accuracy_by_model     = accuracy_by_model,
  recovery_df           = recovery_df,
  n_reps                = opt$n_reps,
  row_totals            = row_totals,
  timestamp             = Sys.time()
), summary_file)
cat(sprintf("Summary saved to: %s\n\n", summary_file))

# ============================================================================
# PLOTS
# ============================================================================

plot_file <- file.path(opt$output_dir, "recovery_analysis.pdf")
pdf(plot_file, width = 12, height = 10)

base_theme <- theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title  = element_text(face = "bold"),
        legend.position = "bottom")
model_factor <- function(x) factor(x, levels = model_names)

# Plot 1: Confusion Matrix
conf_long <- data.frame(
  Generating = model_factor(rep(model_names, each = n_models)),
  Selected   = model_factor(rep(model_names, times = n_models)),
  Proportion = as.vector(t(confusion_prop))
)

p1 <- ggplot(conf_long, aes(x = Selected, y = Generating, fill = Proportion)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", Proportion)), size = 5) +
  scale_fill_gradient2(low = "white", high = "darkblue", mid = "lightblue",
                       midpoint = 0.5, limits = c(0, 1)) +
  labs(title = "Model Recovery: Confusion Matrix",
       subtitle = sprintf("Overall Accuracy: %.1f%% | N = %d simulations",
                          100 * overall_accuracy, sum(row_totals)),
       x = "Selected Model (via LOO)", y = "Generating Model") +
  base_theme
print(p1)

# Plot 2: Inversion Matrix
inv_long <- data.frame(
  Generating  = model_factor(rep(model_names, each = n_models)),
  Selected    = model_factor(rep(model_names, times = n_models)),
  Probability = as.vector(t(inversion_matrix))
)

p2 <- ggplot(inv_long, aes(x = Selected, y = Generating, fill = Probability)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", Probability)), size = 5) +
  scale_fill_gradient2(low = "white", high = "darkgreen", mid = "lightgreen",
                       midpoint = 0.5, limits = c(0, 1)) +
  labs(title = "Inversion Matrix: P(Generated by X | Selected Y)",
       x = "Selected Model", y = "Generating Model") +
  base_theme
print(p2)

# Plot 3: Accuracy by model with CI
p3 <- ggplot(accuracy_by_model,
             aes(x = gen_model_name, y = accuracy)) +
  geom_col(fill = "steelblue", width = 0.6) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * accuracy)), vjust = -0.5, size = 4) +
  scale_y_continuous(limits = c(0, 1.05), labels = scales::percent) +
  labs(title = "Recovery Accuracy by Generating Model",
       subtitle = sprintf("Error bars = 95%% binomial CI | N = %d simulations", nrow(recovery_df)),
       x = "Generating Model", y = "Recovery Accuracy") +
  base_theme
print(p3)

dev.off()
cat(sprintf("Plots saved to: %s\n", plot_file))

# ============================================================================
# LATEX TABLE
# ============================================================================

latex_file <- file.path(opt$output_dir, "recovery_table.tex")
lines <- c(
  "\\begin{table}[h]", "\\centering",
  sprintf("\\caption{Model Recovery Results (N = %d simulations)}", nrow(recovery_df)),
  "\\begin{tabular}{lcccc|cc}", "\\hline",
  "& \\multicolumn{4}{c|}{Selected Model} & & \\\\",
  "Generating Model & Full & No Weights & No Distortion & Null & Recovery\\% & N \\\\",
  "\\hline"
)
for (i in seq_len(n_models)) {
  cells <- paste(sprintf("%.2f", confusion_prop[i, ]), collapse = " & ")
  acc   <- 100 * confusion_matrix[i, i] / max(row_totals[i], 1)
  lines <- c(lines, sprintf("%s & %s & %.1f\\%% & %d \\\\",
                            model_names[i], cells, acc, row_totals[i]))
}
lines <- c(lines, "\\hline",
  sprintf("Overall & \\multicolumn{4}{c|}{} & %.1f\\%% & %d \\\\",
          100 * overall_accuracy, sum(row_totals)),
  "\\hline", "\\end{tabular}", "\\end{table}")
writeLines(lines, latex_file)
cat(sprintf("LaTeX table saved to: %s\n", latex_file))

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("\n=== SUMMARY ===\n")
cat(sprintf("Total compared   : %d\n", nrow(recovery_df)))
cat(sprintf("Overall accuracy : %.1f%%\n", 100 * overall_accuracy))
for (i in seq_len(n_models))
  cat(sprintf("  %-15s %d/%d (%.1f%%)\n", model_names[i],
              confusion_matrix[i, i], row_totals[i],
              100 * confusion_matrix[i, i] / max(row_totals[i], 1)))
cat("\nDone!\n")
