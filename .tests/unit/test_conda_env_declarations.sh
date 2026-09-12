#!/usr/bin/env bash
# Structural audit: every rule that invokes Rscript/python must declare a conda: environment,
# otherwise --use-conda silently falls back to whatever's already on PATH (typically just the
# launcher's own snakemake-env, which has no R/Python packages installed) - confirmed directly
# on a real cluster run: generate_genes_on_bins_obj.smk had no conda: directive at all and
# failed with exit 127 ("command not found") the moment Rscript itself couldn't be found.
# This check is deliberately lightweight (text-based), not a real `conda env create` run - see
# README for why a real --use-conda test isn't part of this project's CI.
# Usage: bash .tests/unit/test_conda_env_declarations.sh

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

# Rules intentionally exempt: run_HiCQC activates a manually pre-created conda env
# ("hic_qc") directly in its own shell block, since it isn't a standard, describable-in-yml
# environment - documented and confirmed intentional, not an oversight.
EXEMPT="workflow/rules/HiC/run_HiCQC.smk"

FAIL=0
MISSING=""
for smk in $(find workflow/rules -name "*.smk"); do
    if grep -qE "Rscript|python3|python -m|python " "$smk" 2>/dev/null; then
        if ! grep -q "conda:" "$smk"; then
            if [ "$smk" != "$EXEMPT" ]; then
                MISSING="$MISSING $smk"
                FAIL=1
            fi
        fi
    fi
done

if [ -n "$MISSING" ]; then
    echo "  [FAIL] rules invoking Rscript/python with no conda: directive:$MISSING"
else
    echo "  [OK] every Rscript/python-invoking rule (except the documented run_HiCQC exception) declares a conda: environment"
fi

echo
if [ "$FAIL" -eq 0 ]; then
    echo "ALL TESTS PASSED"
    exit 0
else
    echo "SOME TESTS FAILED"
    exit 1
fi
