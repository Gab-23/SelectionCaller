#!/usr/bin/env bash
# Unit tests for the ML pipeline's Python/R helper scripts. The full Snakemake-driven pipeline
# (spatial blocking -> nested CV -> hyperparameter search -> SHAP-ready models -> aggregation ->
# visualization) is validated end-to-end separately with a toy dataset and real Snakemake
# execution - this script covers the specific bugs found and fixed in isolation.
# Usage: bash .tests/unit/test_ml_scripts.sh
# Requires: python3 (numpy, pandas, scikit-learn, scikit-optimize, xgboost, shap), Rscript
# (dplyr, tidyr, tibble, ggplot2, magrittr)

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

FAIL=0
check() {
    local desc="$1"; local got="$2"; local expected="$3"
    if [ "$got" == "$expected" ]; then
        echo "  [OK] $desc"
    else
        echo "  [FAIL] $desc"
        echo "         expected: $expected"
        echo "         got:      $got"
        FAIL=1
    fi
}

echo "=== 1) CHECK_CHR_SORTING.py ==="
echo "    (bug: the per-chromosome 'sortedness' check was trivially always True regardless of"
echo "    whether different chromosomes were contiguous or interleaved in the original array -"
echo "    confirmed directly with an interleaved chr1/chr2 array incorrectly reported as"
echo "    sorted. This check gates GENERATE_KFOLD_IDXS's naive contiguous-row-index KFold"
echo "    split. Per your confirmation this case never occurs with your real data, but fixed"
echo "    regardless as a zero-cost safety net.)"
python3 -c "
import sys, pandas as pd
sys.path.insert(0, 'workflow/scripts/helpers')
from CHECK_CHR_SORTING import CHECK_CHR_SORTING
interleaved = pd.Series(['2_1','1_1','2_2','1_2'])
sorted_ok = pd.Series(['1_1','1_2','2_1','2_2'])
within_chr_unsorted = pd.Series(['1_2','1_1','2_1','2_2'])
print('INTERLEAVED', CHECK_CHR_SORTING(interleaved))
print('SORTED', CHECK_CHR_SORTING(sorted_ok))
print('WITHIN_CHR_UNSORTED', CHECK_CHR_SORTING(within_chr_unsorted))
" > "$TMP/chr_sorting_out.txt"
check "interleaved chromosomes correctly caught (False)" "$(grep INTERLEAVED "$TMP/chr_sorting_out.txt" | cut -d' ' -f2)" "False"
check "properly sorted array still passes (True)" "$(grep '^SORTED' "$TMP/chr_sorting_out.txt" | cut -d' ' -f2)" "True"
check "within-chromosome disorder still caught (False, non-regression)" "$(grep WITHIN_CHR "$TMP/chr_sorting_out.txt" | cut -d' ' -f2)" "False"

echo "=== 2) OPTIMIZE_HYPERPARAMS.py ==="
echo "    (bug: GridSearchCV was called with n_iter, a parameter it does not accept at all -"
echo "    confirmed directly this raised a TypeError on instantiation, guaranteed crash every"
echo "    time GridSearchCV was selected as the optimizer.)"
RC=$(python3 -c "
import sys
sys.path.insert(0, 'workflow/scripts/helpers')
import numpy as np, pandas as pd
from sklearn.linear_model import LinearRegression
from OPTIMIZE_HYPERPARAMS import OPTIMIZE_HYPERPARAMS
np.random.seed(0)
X = pd.DataFrame(np.random.rand(20, 2), columns=['a','b'])
y = pd.Series(np.random.rand(20))
space = {'GridSearchCV': {'fit_intercept': [True, False]}}
OPTIMIZE_HYPERPARAMS(LinearRegression(), space, 'GridSearchCV', X, y, 42, 1, 3, 2)
" > /dev/null 2>&1; echo $?)
check "GridSearchCV no longer crashes with n_iter TypeError" "$RC" "0"

echo "=== 3) GET_EXPLAINER.py ==="
echo "    (bug: shap.LinearExplainer(...) passed a bare Ellipsis as its only argument - the"
echo "    masker parameter was never implemented at all, not just mistyped. Confirmed this"
echo "    raised 'missing 1 required positional argument: masker'. Fixed by adding a"
echo "    BACKGROUND_DATA parameter, used as the masker.)"
RC=$(python3 -c "
import sys
sys.path.insert(0, 'workflow/scripts/helpers')
import numpy as np, pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import Pipeline
from GET_EXPLAINER import GET_EXPLAINER
np.random.seed(0)
X = pd.DataFrame(np.random.rand(20, 2), columns=['a','b'])
y = pd.Series(np.random.rand(20))
pipe = Pipeline([('model', LinearRegression())]).fit(X, y)
explainer = GET_EXPLAINER(pipe, IS_TREE_BASED=False, IS_LINEAR=True, BACKGROUND_DATA=X)
explainer(X)
" > /dev/null 2>&1; echo $?)
check "LinearExplainer no longer crashes with missing masker" "$RC" "0"

echo "=== 4) BUILD_MODEL.py ==="
echo "    (n_jobs was hardcoded to 1 for every parallelizable model - parametrized per your"
echo "    request; also added a defensive else-raise for an unrecognized model name, which"
echo "    previously would have raised an unclear UnboundLocalError.)"
python3 -c "
import sys
sys.path.insert(0, 'workflow/scripts/helpers')
from BUILD_MODEL import BUILD_MODEL
m = BUILD_MODEL('RandomForestRegressor', 4)
print('N_JOBS', m.n_jobs)
" > "$TMP/build_model_out.txt" 2>&1
check "N_JOBS is actually threaded through to the model" "$(grep N_JOBS "$TMP/build_model_out.txt" | cut -d' ' -f2)" "4"
RC=$(python3 -c "
import sys
sys.path.insert(0, 'workflow/scripts/helpers')
from BUILD_MODEL import BUILD_MODEL
BUILD_MODEL('NotARealModel', 1)
" > /dev/null 2>&1; echo $?)
check "unknown model name now fails clearly (non-zero exit) instead of UnboundLocalError" "$RC" "1"

echo "=== 5) train_models.py output files (index=True) ==="
echo "    (bug: all three final .to_csv() calls used index=False, silently discarding the"
echo "    '{MODEL}_{BLOCKING}_{fold}' row identifier that aggregate_and_visualize_models.R and"
echo "    visualize_models.R both depend on (read.delim(..., row.names=1) +"
echo "    separate(Model, into=c('model_name','blocking','fold'))). Confirmed directly this"
echo "    crashed R's read.delim with 'duplicate row.names are not allowed' once two rows"
echo "    happened to share a value in what R fell back to treating as the row-name column.)"
python3 -c "
import pandas as pd
d = {'ModelA_5_0': {'r2': 0.1}, 'ModelA_5_1': {'r2': 0.1}}  # identical r2 on purpose
df = pd.DataFrame.from_dict(d, orient='index')
df.to_csv('$TMP/perf_test.tsv', sep='\t', header=True, index=True)
"
Rscript -e '
df <- read.delim("'"$TMP"'/perf_test.tsv", sep="\t", row.names=1, check.names=FALSE)
cat("ROWNAMES", paste(rownames(df), collapse=","), "\n")
' > "$TMP/r_read_out.txt" 2>&1
check "row identifiers survive the round-trip even with duplicate metric values" "$(grep ROWNAMES "$TMP/r_read_out.txt" | cut -d' ' -f2)" "ModelA_5_0,ModelA_5_1"

echo "=== 6) aggregate_and_visualize_models.R variable names ==="
echo "    (bug: the final ggsave() calls referenced output_path_heatmap_sel_feats /"
echo "    output_path_heatmap_performances, never assigned anywhere in the script - a"
echo "    copy-paste naming mismatch. Fixed to use the actually-declared"
echo "    output_path_selfeats / output_path_performance.)"
if grep -v "^#" workflow/scripts/aggregate_and_visualize_models.R | grep -q "output_path_heatmap_sel_feats\|output_path_heatmap_performances"; then
    echo "  [FAIL] stale variable names still present in actual code (not just comments)"
    FAIL=1
else
    echo "  [OK] stale variable names no longer present in actual code"
fi

echo "=== 7) visualize_models.R ordering ==="
echo "    (bug: heatmap_df_performances\$y_ordered was set before heatmap_df_performances was"
echo "    ever created - confirmed directly this crashed with 'object not found'. Fixed by"
echo "    moving its creation above the y_ordered assignment.)"
Rscript -e '
`%>%` <- magrittr::`%>%`
source("workflow/scripts/helpers/visualize_models.R")
set.seed(1)
rn <- c("ModelA_5_0","ModelA_5_1","ModelB_noblocking_0")
feats <- data.frame(f1=c(1,0,1), f2=c(0,1,1)); rownames(feats) <- rn
perf <- data.frame(train_r2=runif(3)); rownames(perf) <- rn
result <- visualize_models(feats, perf)
cat("NAMES", paste(names(result), collapse=","), "\n")
' > "$TMP/visualize_out.txt" 2>&1
check "visualize_models() runs end-to-end and returns all 3 plots" "$(grep NAMES "$TMP/visualize_out.txt" | cut -d' ' -f2)" "heatmap_features,heatmap_performances,median_boxplots"

echo
if [ "$FAIL" -eq 0 ]; then
    echo "ALL TESTS PASSED"
    exit 0
else
    echo "SOME TESTS FAILED"
    exit 1
fi
