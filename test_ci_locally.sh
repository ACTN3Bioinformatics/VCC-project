#!/bin/bash
#
# Test CI Locally
# This script simulates the CI tests for the VCC-project locally.
# It checks for critical files, Python imports, code quality, and Snakemake rules.

set -e

# Colors
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

print_success() { echo -e "${GREEN}✓ $1${NC}"; }
print_error() { echo -e "${RED}✗ $1${NC}"; }
print_info() { echo -e "${YELLOW}➜ $1${NC}"; }
print_section() { echo -e "${BLUE}==== $1 ====${NC}"; }

FAILED=0

echo "=========================================="
echo "CI TESTS - LOCAL SIMULATION"
echo "=========================================="
echo ""

# Check if in correct directory
if [ ! -f "Snakefile" ]; then
    print_error "Not in VCC-project root directory"
    exit 1
fi

# Activate environment
print_section "Activating environment"
if mamba env list | grep -q "vcc2025"; then
    source $MAMBA_ROOT_PREFIX/etc/profile.d/mamba.sh
    mamba activate vcc2025
    print_success "Environment activated"
else
    print_error "Environment 'vcc2025' not found"
    exit 1
fi

# Test 1: Structure Check
print_section "Test 1: Project Structure"

CRITICAL_FILES=(
    "Snakefile"
    "environment.yml"
    "config/config.yaml"
    "config/datasets.yaml"
    "workflows/rules/qc.smk"
    "workflows/rules/normalize.smk"
    "workflows/rules/balance.smk"
    "workflows/rules/integrate.smk"
    "workflows/rules/split.smk"
    "workflows/rules/features.smk"
    "workflows/rules/download.smk"
    "scripts/filter_normalize.py"
    "scripts/balance.py"
    "scripts/integration.py"
    "scripts/split_data.py"
    "scripts/feature_engineering.py"
    "scripts/utils.py"
)

for file in "${CRITICAL_FILES[@]}"; do
    if [ -f "$file" ]; then
        print_success "$file"
    else
        print_error "$file MISSING"
        FAILED=1
    fi
done

# Check for deprecated files
if [ -f "workflows/Snakefile" ]; then
    print_error "⚠ workflows/Snakefile exists (should be deprecated)"
    echo "  Run: mv workflows/Snakefile workflows/Snakefile.deprecated"
fi

# Test 2: Python Imports
print_section "Test 2: Python Imports"

python -c "import sys; sys.path.append('scripts'); import utils" 2>/dev/null
if [ $? -eq 0 ]; then
    print_success "utils.py imports OK"
else
    print_error "utils.py import failed"
    FAILED=1
fi

python -c "import scanpy; print(f'scanpy {scanpy.__version__}')" 2>/dev/null
if [ $? -eq 0 ]; then
    print_success "scanpy OK"
else
    print_error "scanpy import failed"
    FAILED=1
fi

python -c "import pandas; import numpy; import snakemake" 2>/dev/null
if [ $? -eq 0 ]; then
    print_success "Core packages OK"
else
    print_error "Some packages missing"
    FAILED=1
fi

# Test 3: Flake8 Syntax Check
print_section "Test 3: Code Quality (flake8)"

if command -v flake8 &> /dev/null; then
    flake8 scripts/ --count --select=E9,F63,F7,F82 --show-source --statistics
    if [ $? -eq 0 ]; then
        print_success "No syntax errors"
    else
        print_error "Syntax errors found"
        FAILED=1
    fi
else
    print_info "flake8 not installed, skipping"
fi

# Test 4: Snakemake Test Rule
print_section "Test 4: Snakemake Test Rule"

snakemake test --cores 1 2>&1 | tee /tmp/snakemake_test.log
if [ $? -eq 0 ]; then
    print_success "Snakemake test rule OK"
else
    print_error "Snakemake test rule failed"
    echo "See: /tmp/snakemake_test.log"
    FAILED=1
fi

# Test 5: Snakemake Dry-run
print_section "Test 5: Snakemake Dry-run"

snakemake -n --cores 1 2>&1 | tee /tmp/snakemake_dryrun.log
if [ $? -eq 0 ]; then
    print_success "Snakemake dry-run OK"
else
    print_error "Snakemake dry-run failed"
    echo "See: /tmp/snakemake_dryrun.log"
    FAILED=1
fi

# Test 6: DAG Construction
print_section "Test 6: DAG Construction"

snakemake --dag > /tmp/dag.dot 2>&1
if [ $? -eq 0 ]; then
    print_success "DAG construction OK"
    echo "  DAG saved to: /tmp/dag.dot"
else
    print_error "DAG construction failed"
    FAILED=1
fi

# Summary
echo ""
echo "=========================================="
print_section "SUMMARY"
echo "=========================================="

if [ $FAILED -eq 0 ]; then
    print_success "All tests PASSED ✨"
    echo ""
    echo "Your code is ready for CI/CD!"
    echo ""
    echo "Next steps:"
    echo "  1. git add ."
    echo "  2. git commit -m 'Your message'"
    echo "  3. git push origin main"
    echo ""
    exit 0
else
    print_error "Some tests FAILED"
    echo ""
    echo "Please fix the issues above before pushing."
    echo ""
    exit 1
fi
