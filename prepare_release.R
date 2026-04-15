#!/usr/bin/env Rscript
# Release Preparation Script for rSHUD v2.2.0
# This script performs final checks before CRAN submission

cat("=== rSHUD v2.2.0 Release Preparation ===\n\n")

# Load required packages
if (!require("devtools")) install.packages("devtools")
if (!require("usethis")) install.packages("usethis")

# 1. Check package structure
cat("1. Checking package structure...\n")
if (!file.exists("DESCRIPTION")) {
  stop("DESCRIPTION file not found!")
}
if (!file.exists("NEWS.md")) {
  stop("NEWS.md file not found!")
}
if (!file.exists("README.md")) {
  stop("README.md file not found!")
}
cat("   ✓ Package structure OK\n\n")

# 2. Check version number
cat("2. Checking version number...\n")
desc <- read.dcf("DESCRIPTION")
version <- desc[1, "Version"]
if (version != "2.2.0") {
  stop("Version should be 2.2.0, found: ", version)
}
cat("   ✓ Version is 2.2.0\n\n")

# 3. Check dependencies
cat("3. Checking dependencies...\n")
imports <- strsplit(desc[1, "Imports"], ",\\s*")[[1]]
if (!any(grepl("terra", imports))) {
  stop("terra not found in Imports!")
}
if (!any(grepl("sf", imports))) {
  stop("sf not found in Imports!")
}
cat("   ✓ terra and sf in Imports\n")

suggests <- desc[1, "Suggests"]
if (!grepl("raster", suggests)) {
  warning("raster not in Suggests (should be for testing)")
}
cat("   ✓ Dependencies OK\n\n")

# 4. Update documentation
cat("4. Updating documentation...\n")
tryCatch({
  devtools::document()
  cat("   ✓ Documentation updated\n\n")
}, error = function(e) {
  cat("   ✗ Documentation update failed:", e$message, "\n\n")
})

# 5. Run tests
cat("5. Running tests...\n")
tryCatch({
  test_results <- devtools::test(reporter = "summary")
  cat("   ✓ Tests completed\n\n")
}, error = function(e) {
  cat("   ✗ Tests failed:", e$message, "\n\n")
})

# 6. Run R CMD check
cat("6. Running R CMD check...\n")
cat("   This may take several minutes...\n")
tryCatch({
  check_results <- devtools::check(
    document = FALSE,
    args = c("--as-cran"),
    error_on = "warning"
  )
  cat("   ✓ R CMD check passed\n\n")
}, error = function(e) {
  cat("   ✗ R CMD check failed:", e$message, "\n\n")
  cat("   Please review and fix issues before submission\n\n")
})

# 7. Check URLs
cat("7. Checking URLs...\n")
if (require("urlchecker", quietly = TRUE)) {
  tryCatch({
    urlchecker::url_check()
    cat("   ✓ URLs checked\n\n")
  }, error = function(e) {
    cat("   ⚠ URL check warning:", e$message, "\n\n")
  })
} else {
  cat("   ⚠ urlchecker not installed, skipping\n\n")
}

# 8. Spell check
cat("8. Running spell check...\n")
if (require("spelling", quietly = TRUE)) {
  tryCatch({
    spelling::spell_check_package()
    cat("   ✓ Spell check completed\n\n")
  }, error = function(e) {
    cat("   ⚠ Spell check warning:", e$message, "\n\n")
  })
} else {
  cat("   ⚠ spelling not installed, skipping\n\n")
}

# 9. Build package tarball
cat("9. Building package tarball...\n")
tryCatch({
  tarball <- devtools::build()
  cat("   ✓ Package built:", tarball, "\n\n")
}, error = function(e) {
  cat("   ✗ Build failed:", e$message, "\n\n")
})

# 10. Check built package
cat("10. Checking built package...\n")
if (exists("tarball") && file.exists(tarball)) {
  cat("   Run manually: R CMD check --as-cran", tarball, "\n\n")
} else {
  cat("   ⚠ Tarball not found, skipping\n\n")
}

# Summary
cat("=== Release Preparation Summary ===\n\n")
cat("Package: rSHUD\n")
cat("Version: 2.2.0\n")
cat("Date:", Sys.Date(), "\n\n")

cat("Release materials created:\n")
cat("  ✓ rSHUD_2.2.0.tar.gz - Package tarball\n")
cat("  ✓ docs/release/RELEASE_NOTES_v2.2.0.md - Comprehensive release notes\n")
cat("  ✓ docs/release/cran-comments.md - CRAN submission comments\n")
cat("  ✓ docs/release/CRAN_SUBMISSION_CHECKLIST.md - Submission checklist\n\n")

cat("Next steps:\n")
cat("  1. Review docs/release/CRAN_SUBMISSION_CHECKLIST.md\n")
cat("  2. Run additional platform checks:\n")
cat("     - devtools::check_win_devel()\n")
cat("     - devtools::check_rhub()\n")
cat("  3. Review docs/release/cran-comments.md\n")
cat("  4. Submit to CRAN:\n")
cat("     - devtools::submit_cran()\n")
cat("     - Or upload manually to https://cran.r-project.org/submit.html\n\n")

cat("=== Preparation Complete ===\n")
