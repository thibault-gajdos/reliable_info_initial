#!/usr/bin/env Rscript

# ============================================================================
# CLEAN renv SETUP FOR CLUSTER (SOURCE INSTALL, USER LIBRARY)
# ============================================================================

cat("=== Starting clean renv setup ===\n\n")

# ---- 1) FIX ENVIRONMENT (CRITICAL) ----

# Use user-owned library (NOT /usr/local)
Sys.setenv(RENV_PATHS_LIBRARY = "~/renv/library")

# Disable staging (fixes cross-device /tmp issues on clusters)
Sys.setenv(RENV_CONFIG_INSTALL_STAGED = "FALSE")

# Force source installs (Linux clusters don't support binaries)
options(pkgType = "source")

# Avoid aggressive parallel compilation (safer on clusters)
Sys.setenv(MAKEFLAGS = "-j1")

# Optional: ensure reproducibility
options(repos = c(CRAN = "https://cloud.r-project.org"))

# ---- 2) LOAD renv ----

suppressMessages(library(renv))

cat("Library paths after setup:\n")
print(.libPaths())
cat("\n")

# ---- 3) CLEAN PREVIOUS BROKEN STATE ----

cat("=== Cleaning previous renv (if any) ===\n")

tryCatch({
  renv::deactivate()
}, error = function(e) {})

if (dir.exists("renv")) {
  unlink("renv", recursive = TRUE, force = TRUE)
  cat("Removed existing renv folder\n")
}

if (file.exists("renv.lock")) {
  unlink("renv.lock")
  cat("Removed existing renv.lock\n")
}

cat("\n")

# ---- 4) INITIALIZE renv ----

cat("=== Initializing renv ===\n")
renv::init(bare = TRUE)

cat("renv initialized\n\n")

# ---- 5) INSTALL CORE TOOLCHAIN PACKAGES FIRST ----

cat("=== Installing core packages ===\n")

core_pkgs <- c("rlang", "Rcpp", "RcppEigen")

for (pkg in core_pkgs) {
  cat(sprintf("\nInstalling %s ...\n", pkg))
  
  tryCatch({
    renv::install(pkg)
    
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("%s installed but cannot be loaded", pkg))
    }
    
    cat(sprintf("✓ %s OK\n", pkg))
    
  }, error = function(e) {
    cat(sprintf("✗ FAILED: %s\n", e$message))
    quit(status = 1)
  })
}

cat("\nCore packages installed successfully\n\n")

# ---- 6) INSTALL REMAINING PACKAGES ----

cat("=== Installing remaining packages ===\n")

pkgs <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "optparse",
  "loo",
  "stan-dev/cmdstanr"
)

for (pkg in pkgs) {
  cat(sprintf("\nInstalling %s ...\n", pkg))
  
  tryCatch({
    renv::install(pkg)
    
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("%s installed but cannot be loaded", pkg))
    }
    
    cat(sprintf("✓ %s OK\n", pkg))
    
  }, error = function(e) {
    cat(sprintf("✗ FAILED: %s\n", pkg))
    cat(sprintf("Reason: %s\n", e$message))
  })
}

# ---- 7) SNAPSHOT ----

cat("\n=== Creating renv snapshot ===\n")

tryCatch({
  renv::snapshot()
  cat("✓ Snapshot created\n")
}, error = function(e) {
  cat(sprintf("Snapshot failed: %s\n", e$message))
})

# ---- 8) FINAL CHECK ----

cat("\n=== Final package check ===\n")

all_packages <- c("rlang", "Rcpp", "RcppEigen",
                  "dplyr", "tidyr", "ggplot2",
                  "optparse", "loo", "rstan")

all_ok <- TRUE

for (pkg in all_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("✓ %-12s %s\n", pkg, as.character(packageVersion(pkg))))
  } else {
    cat(sprintf("✗ %-12s NOT INSTALLED\n", pkg))
    all_ok <- FALSE
  }
}

# ---- 9) RESULT ----

if (all_ok) {
  cat("\n=== SUCCESS ===\n")
  cat("All packages installed correctly.\n")
} else {
  cat("\n=== PARTIAL FAILURE ===\n")
  cat("Some packages failed.\n")
  cat("If rstan failed, likely cause is missing compiler toolchain.\n")
}

cat("\n=== Done ===\n")
