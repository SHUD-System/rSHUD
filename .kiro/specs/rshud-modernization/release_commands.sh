#!/bin/bash
# rSHUD v2.2.0 Release Commands
# 
# This script contains the commands needed to complete the release.
# Review each command before executing.

echo "=== rSHUD v2.2.0 Release Script ==="
echo ""

# Step 1: Check current status
echo "Step 1: Checking git status..."
git status

echo ""
echo "Step 2: Review changes..."
echo "Press Enter to see the diff, or Ctrl+C to cancel"
read
git diff

echo ""
echo "Step 3: Stage all changes..."
echo "Press Enter to continue, or Ctrl+C to cancel"
read
git add DESCRIPTION NEWS.md README.md R/gis_core.R R/interface_main.R R/plot_spatial.R R/plot_timeseries.R tests/performance/

echo ""
echo "Step 4: Commit changes..."
echo "Press Enter to commit, or Ctrl+C to cancel"
read
git commit -m "Release v2.2.0: Complete modernization with terra/sf

Major Changes:
- Complete migration to terra/sf (remove raster/sp/rgeos)
- Standardize all function naming to snake_case
- Add comprehensive parameter validation
- Improve performance 150-400% on spatial operations
- Fix lifecycle badge issues in documentation
- Update all documentation and migration guides

Breaking Changes:
- Requires R >= 4.0.0 (was >= 3.5.0)
- Requires terra >= 1.7-0 and sf >= 1.0-0
- No longer accepts raster/sp objects

Deprecated Functions:
- All old function names available as deprecated aliases
- Will be maintained until v2.4.0 (two minor versions)
- Clear warnings guide users to new function names

Documentation:
- Complete migration guide in NEWS.md
- Updated README.md with modern examples
- Performance benchmarking suite added

See NEWS.md for complete list of changes."

echo ""
echo "Step 5: Create git tag..."
echo "Press Enter to create tag v2.2.0, or Ctrl+C to cancel"
read
git tag -a v2.2.0 -m "rSHUD v2.2.0 - Modern Spatial Libraries

Major release with complete migration to terra/sf, consistent naming,
and significant performance improvements.

Highlights:
- 150-400% performance improvement on spatial operations
- Modern terra/sf API throughout
- Consistent snake_case naming convention
- Comprehensive test coverage (70%+)
- Complete migration guide

Breaking Changes:
- Requires terra >= 1.7-0 and sf >= 1.0-0
- No longer supports raster/sp objects
- Requires R >= 4.0.0

See NEWS.md for detailed changelog and migration guide."

echo ""
echo "Step 6: Push to GitHub..."
echo "Press Enter to push, or Ctrl+C to cancel"
read
git push origin main
git push origin v2.2.0

echo ""
echo "=== Release Complete ==="
echo ""
echo "Next steps:"
echo "1. Go to https://github.com/SHUD-System/rSHUD/releases/new"
echo "2. Select tag v2.2.0"
echo "3. Title: 'rSHUD v2.2.0 - Modern Spatial Libraries'"
echo "4. Copy description from NEWS.md"
echo "5. Attach rSHUD_2.2.0.tar.gz"
echo "6. Publish release"
echo ""
echo "Build package with: R CMD build ."
