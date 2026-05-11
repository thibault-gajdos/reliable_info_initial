#!/usr/bin/env Rscript
# ============================================================================
# Plot Parameter Recovery — Choice-Only Model
# ============================================================================

rm(list = ls(all = TRUE))
setwd("~/reliable_info/param_recovery_choice")
renv::load()

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

results_dir <- '../results/choice'
plot_dir    <- '../figures/choice'
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# 1. LOAD ALL RESULTS
# ============================================================================

files <- list.files(results_dir, pattern = "^recover_\\d+\\.rds$", full.names = TRUE)
cat(sprintf("Found %d result files.\n", length(files)))
if (length(files) == 0) stop("No result files found.")

param_info <- list(
  alpha = list(sim = "Simulated a",    fit = "Fitted a"),
  beta  = list(sim = "Simulated b",    fit = "Fitted b"),
  w1    = list(sim = "Simulated w(1)", fit = "Fitted w(1)"),
  w2    = list(sim = "Simulated w(2)", fit = "Fitted w(2)"),
  w3    = list(sim = "Simulated w(3)", fit = "Fitted w(3)"),
  w4    = list(sim = "Simulated w(4)", fit = "Fitted w(4)"),
  w5    = list(sim = "Simulated w(5)", fit = "Fitted w(5)")
)
param_order   <- names(param_info)
panel_letters <- LETTERS[seq_along(param_order)]

group_rows <- list()
indiv_rows <- list()
diag_rows  <- list()

for (f in files) {
  res <- readRDS(f)

  ## Group-level
  gf <- res$group_fitted
  group_rows[[f]] <- tibble(
    k         = res$k,
    param     = res$param_names,
    simulated = res$group_sim,
    estimated = gf$mean,
    q5        = gf$q5,
    q95       = gf$q95
  )

  ## Individual-level
  ## indiv_sim is a flat vector: alpha[1..N], beta[1..N], w[1,1]..w[1,5], w[2,1]..w[2,5], ...
  ## indiv_fitted follows the same order
  inf <- res$indiv_fitted
  N   <- length(res$indiv_sim) / 7  # 7 params per subject

  param_vec <- c(
    rep("alpha", N),
    rep("beta", N),
    rep(c("w1","w2","w3","w4","w5"), times = N)
  )
  subject_vec <- c(
    seq_len(N),
    seq_len(N),
    rep(seq_len(N), each = 5)
  )

  indiv_rows[[f]] <- tibble(
    k         = res$k,
    subject   = subject_vec,
    param     = param_vec,
    simulated = res$indiv_sim,
    estimated = inf$mean
  )

  ## Diagnostics
  diag_rows[[f]] <- tibble(
    k           = res$k,
    n_divergent = res$n_divergent,
    n_max_td    = res$n_max_td,
    max_rhat    = res$max_rhat,
    min_ess     = res$min_ess
  )
}

df_group <- bind_rows(group_rows)
df_indiv <- bind_rows(indiv_rows)
df_diag  <- bind_rows(diag_rows)

## Filter problematic fits
good_fits <- df_diag %>% filter(max_rhat < 1.1, n_divergent < 50)
cat(sprintf("Keeping %d / %d fits (Rhat < 1.1 & divergences < 50)\n",
            nrow(good_fits), nrow(df_diag)))
df_group <- df_group %>% filter(k %in% good_fits$k)
df_indiv <- df_indiv %>% filter(k %in% good_fits$k)

# ============================================================================
# 2. DIAGNOSTICS SUMMARY
# ============================================================================

cat("\n=== Diagnostics Summary ===\n")
cat(sprintf("  Fits with divergences:     %d / %d\n",
            sum(df_diag$n_divergent > 0), nrow(df_diag)))
cat(sprintf("  Fits with Rhat > 1.05:     %d / %d\n",
            sum(df_diag$max_rhat > 1.05), nrow(df_diag)))
cat(sprintf("  Fits with ESS < 400:       %d / %d\n",
            sum(df_diag$min_ess < 400), nrow(df_diag)))
cat("\n")

# ============================================================================
# HELPER: format p-value
# ============================================================================

format_p <- function(p) {
  if (p < 2.2e-16) return("italic(p) < 2.2e-16")
  sprintf("italic(p) == %.1e", p)
}

# ============================================================================
# HELPER: make one panel (matching reference style)
# ============================================================================

make_panel <- function(df, par, letter, show_errorbars = FALSE, axis_lim = NULL) {
  dd <- df %>% filter(param == par)

  if (is.null(axis_lim)) {
    all_vals <- c(dd$simulated, dd$estimated)
    if (show_errorbars && all(c("q5", "q95") %in% names(dd)))
      all_vals <- c(all_vals, dd$q5, dd$q95)
    axis_lim <- range(all_vals, na.rm = TRUE)
    axis_lim <- axis_lim + c(-0.08, 0.08) * diff(axis_lim)
  }

  ct <- cor.test(dd$simulated, dd$estimated)
  label_text <- sprintf("italic(R) == %.2f*','~%s", ct$estimate, format_p(ct$p.value))

  p <- ggplot(dd, aes(x = simulated, y = estimated)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "black") +
    geom_smooth(method = "lm", se = FALSE, colour = "red", linewidth = 0.8)

  if (show_errorbars && all(c("q5", "q95") %in% names(dd))) {
    p <- p + geom_linerange(aes(ymin = q5, ymax = q95), linewidth = 0.3, colour = "black")
  }

  p <- p +
    geom_point(size = ifelse(show_errorbars, 1.5, 0.5),
               alpha = ifelse(show_errorbars, 0.8, 0.3),
               colour = "black") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.08, vjust = 1.5,
             label = label_text, parse = TRUE, size = 4.5) +
    scale_x_continuous(limits = axis_lim, expand = c(0, 0)) +
    scale_y_continuous(limits = axis_lim, expand = c(0, 0)) +
    coord_fixed() +
    labs(x = param_info[[par]]$sim,
         y = param_info[[par]]$fit,
         tag = letter) +
    theme_classic(base_size = 14) +
    theme(
      plot.tag = element_text(face = "bold", size = 13),
      axis.title = element_text(size = 13),
      plot.margin = margin(5, 10, 5, 5)
    )

  return(p)
}

compute_shared_lim <- function(df, params, show_errorbars = FALSE) {
  dd <- df %>% filter(param %in% params)
  all_vals <- c(dd$simulated, dd$estimated)
  if (show_errorbars && all(c("q5", "q95") %in% names(dd)))
    all_vals <- c(all_vals, dd$q5, dd$q95)
  lim <- range(all_vals, na.rm = TRUE)
  lim + c(-0.05, 0.05) * diff(lim)
}

# ============================================================================
# 3. GROUP-LEVEL RECOVERY PLOT
# ============================================================================

cat("=== Group-Level Recovery ===\n")
group_stats <- df_group %>%
  group_by(param) %>%
  summarise(r = cor(simulated, estimated),
            bias = mean(estimated - simulated),
            coverage = mean(simulated >= q5 & simulated <= q95),
            .groups = "drop")
print(as.data.frame(group_stats), digits = 3)
cat("\n")

w_params <- c("w1","w2","w3","w4","w5")
w_lim_group <- compute_shared_lim(df_group, w_params, show_errorbars = TRUE)

panels_group <- lapply(seq_along(param_order), function(i) {
  par <- param_order[i]
  lim <- if (par %in% w_params) w_lim_group else NULL
  make_panel(df_group, par, panel_letters[i], show_errorbars = TRUE, axis_lim = lim)
})

## Layout: 2 rows — top: alpha, beta, w1, w2  |  bottom: w3, w4, w5
p_group <- (panels_group[[1]] | panels_group[[2]] | panels_group[[3]] | panels_group[[4]]) /
           (panels_group[[5]] | panels_group[[6]] | panels_group[[7]] | plot_spacer())

ggsave(file.path(plot_dir, "recovery_group.pdf"), p_group, width = 14, height = 7)
ggsave(file.path(plot_dir, "recovery_group.png"), p_group, width = 14, height = 7, dpi = 300)
cat("Group plot saved.\n")

# ============================================================================
# 4. INDIVIDUAL-LEVEL RECOVERY PLOT
# ============================================================================

cat("=== Individual-Level Recovery ===\n")
indiv_stats <- df_indiv %>%
  group_by(param) %>%
  summarise(r = cor(simulated, estimated),
            bias = mean(estimated - simulated),
            .groups = "drop")
print(as.data.frame(indiv_stats), digits = 3)
cat("\n")

w_lim_indiv <- compute_shared_lim(df_indiv, w_params, show_errorbars = FALSE)

panels_indiv <- lapply(seq_along(param_order), function(i) {
  par <- param_order[i]
  lim <- if (par %in% w_params) w_lim_indiv else NULL
  make_panel(df_indiv, par, panel_letters[i], show_errorbars = FALSE, axis_lim = lim)
})

p_indiv <- (panels_indiv[[1]] | panels_indiv[[2]] | panels_indiv[[3]] | panels_indiv[[4]]) /
           (panels_indiv[[5]] | panels_indiv[[6]] | panels_indiv[[7]] | plot_spacer())

ggsave(file.path(plot_dir, "recovery_individual.pdf"), p_indiv, width = 14, height = 7)
ggsave(file.path(plot_dir, "recovery_individual.png"), p_indiv, width = 14, height = 7, dpi = 300)
cat("Individual plot saved.\n")

# ============================================================================
# 5. DIAGNOSTICS PLOT
# ============================================================================

p_diag <- (
  ggplot(df_diag, aes(x = max_rhat)) +
    geom_histogram(bins = 50, fill = "grey30") +
    geom_vline(xintercept = 1.05, linetype = "dashed", colour = "red") +
    labs(x = "Max Rhat", y = "Count", tag = "A") +
    theme_classic(base_size = 11) +
    theme(plot.tag = element_text(face = "bold"))
) + (
  ggplot(df_diag, aes(x = min_ess)) +
    geom_histogram(bins = 50, fill = "grey30") +
    geom_vline(xintercept = 400, linetype = "dashed", colour = "red") +
    labs(x = "Min bulk ESS", y = "Count", tag = "B") +
    theme_classic(base_size = 11) +
    theme(plot.tag = element_text(face = "bold"))
) + (
  ggplot(df_diag, aes(x = n_divergent)) +
    geom_histogram(bins = 50, fill = "grey30") +
    labs(x = "N divergent", y = "Count", tag = "C") +
    theme_classic(base_size = 11) +
    theme(plot.tag = element_text(face = "bold"))
)

ggsave(file.path(plot_dir, "recovery_diagnostics.pdf"), p_diag,
       width = 14, height = 4)
ggsave(file.path(plot_dir, "recovery_diagnostics.png"), p_diag,
       width = 14, height = 4, dpi = 300)
cat("Diagnostics plot saved.\n")

# ============================================================================
# 6. PARAMETER CORRELATION PLOT (within-fit correlations)
# ============================================================================

cat("Building parameter correlation plot...\n")

param_cols <- c("alpha","beta","w1","w2","w3","w4","w5")
n_params <- length(param_cols)

df_indiv_wide <- df_indiv %>%
  select(k, subject, param, estimated) %>%
  pivot_wider(names_from = param, values_from = estimated)

all_k <- unique(df_indiv_wide$k)
cor_mats <- array(NA_real_, dim = c(n_params, n_params, length(all_k)))

for (i in seq_along(all_k)) {
  dd <- df_indiv_wide %>% filter(k == all_k[i])
  mat <- as.matrix(dd[, param_cols])
  cor_mats[, , i] <- cor(mat, use = "pairwise.complete.obs")
}

mean_cor <- apply(cor_mats, c(1, 2), mean, na.rm = TRUE)
rownames(mean_cor) <- param_cols
colnames(mean_cor) <- param_cols

display_labels <- c("a","b","w(1)","w(2)","w(3)","w(4)","w(5)")

cor_long <- expand.grid(
  x = param_cols,
  y = param_cols,
  stringsAsFactors = FALSE
)
cor_long$r <- as.vector(mean_cor)
cor_long$x <- factor(cor_long$x, levels = param_cols)
cor_long$y <- factor(cor_long$y, levels = rev(param_cols))

p_corr <- ggplot(cor_long, aes(x = x, y = y, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", r)), size = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1),
                       name = "Mean r") +
  scale_x_discrete(labels = display_labels) +
  scale_y_discrete(labels = rev(display_labels)) +
  labs(title = sprintf("Mean parameter correlations (across %d fits)", length(all_k))) +
  coord_fixed() +
  theme_minimal(base_size = 11) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

ggsave(file.path(plot_dir, "param_correlations.pdf"), p_corr,
       width = 7, height = 6)
ggsave(file.path(plot_dir, "param_correlations.png"), p_corr,
       width = 7, height = 6, dpi = 300)

cat("\nMean within-fit correlation matrix:\n")
print(round(mean_cor, 3))
cat("\n")

cat("Parameter correlation plot saved.\n")

cat(sprintf("\nAll figures saved to: %s\n", plot_dir))
