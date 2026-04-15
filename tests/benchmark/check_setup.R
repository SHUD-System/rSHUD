#!/usr/bin/env Rscript
# Check if benchmark environment is properly set up

cat("\n=== Benchmark Setup Check ===\n\n")

# Check required packages
required_pkgs <- c("microbenchmark", "terra", "sf")
optional_pkgs <- c("raster", "sp", "rgeos")

cat("Checking required packages:\n")
for (pkg in required_pkgs) {
  installed <- requireNamespace(pkg, quietly = TRUE)
  status <- if (installed) "✓ Installed" else "✗ Missing"
  cat(sprintf("  %-20s: %s\n", pkg, status))
  
  if (!installed) {
    cat(sprintf("    Install with: install.packages('%s')\n", pkg))
  }
}

cat("\nChecking optional packages (for legacy comparison):\n")
for (pkg in optional_pkgs) {
  installed <- requireNamespace(pkg, quietly = TRUE)
  status <- if (installed) "✓ Installed" else "○ Not installed"
  cat(sprintf("  %-20s: %s\n", pkg, status))
}

# Check if all required packages are available
all_required <- all(sapply(required_pkgs, requireNamespace, quietly = TRUE))
any_optional <- any(sapply(optional_pkgs, requireNamespace, quietly = TRUE))

cat("\n")
if (all_required) {
  cat("✓ All required packages are installed.\n")
  cat("  You can run benchmarks with: source('tests/benchmark/bench-raster.R')\n")
  
  if (any_optional) {
    cat("\n✓ Legacy packages available - full comparison benchmarks can run.\n")
  } else {
    cat("\n○ Legacy packages not available - benchmarks will test modern implementations only.\n")
    cat("  To enable comparison benchmarks, install:\n")
    cat("    install.packages(c('raster', 'sp', 'rgeos'))\n")
  }
} else {
  cat("✗ Some required packages are missing.\n")
  cat("  Install missing packages before running benchmarks.\n")
}

# Check R version
cat("\nR version:\n")
cat(sprintf("  %s\n", R.version.string))

# Check system info
cat("\nSystem information:\n")
cat(sprintf("  OS: %s\n", Sys.info()["sysname"]))
cat(sprintf("  Platform: %s\n", R.version$platform))

cat("\n=== Setup Check Complete ===\n\n")
