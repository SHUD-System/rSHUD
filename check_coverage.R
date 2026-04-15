#!/usr/bin/env Rscript

# Check test coverage for rSHUD package
cat("Checking test coverage for rSHUD...\n\n")

# Check if covr is installed
if (!requireNamespace("covr", quietly = TRUE)) {
  cat("Installing covr package...\n")
  install.packages("covr", repos = "https://cloud.r-project.org")
}

# Load covr
library(covr)

# Calculate package coverage
cat("Calculating coverage (this may take a few minutes)...\n")
cov <- package_coverage(type = "tests", quiet = FALSE)

# Print summary
cat("\n=== Coverage Summary ===\n")
print(cov)

# Get coverage percentage
pct <- percent_coverage(cov)
cat("\nOverall Coverage:", round(pct, 2), "%\n")

# Check if meets target
target <- 70
if (pct >= target) {
  cat("✓ Coverage meets target of", target, "%\n")
} else {
  cat("✗ Coverage below target of", target, "%\n")
  cat("  Need to improve by:", round(target - pct, 2), "%\n")
}

# Show coverage by file
cat("\n=== Coverage by File ===\n")
file_cov <- coverage_to_list(cov)
for (file in names(file_cov)) {
  file_pct <- file_cov[[file]]$coverage
  cat(sprintf("%-40s %6.2f%%\n", basename(file), file_pct))
}

# Identify files with low coverage
cat("\n=== Files with Coverage < 50% ===\n")
low_cov_files <- names(file_cov)[sapply(file_cov, function(x) x$coverage < 50)]
if (length(low_cov_files) > 0) {
  for (file in low_cov_files) {
    cat("  -", basename(file), "\n")
  }
} else {
  cat("  None - all files have >= 50% coverage\n")
}

cat("\nDone!\n")
