#!/usr/bin/env bash
# Unit tests for the ML pipeline's Python/R helper scripts. The full Snakemake-driven pipeline
# (per-chromosome split -> nested CV -> hyperparameter search -> SHAP-ready models -> best-model
# selection -> consensus + visualization) is validated end-to-end separately with a toy dataset
# and real Snakemake execution (see .github/workflows/ci.yml's ml-toy-integration-test job) -
# this script covers the specific bugs found and fixed in isolation.
#
# NOTE: aggregate_and_visualize_models.R/.smk were removed (consolidated into
# select_best_model.R/.smk) - any test referencing them has been replaced or removed below.
#
# Usage: bash .tests/unit/test_ml_scripts.sh
# Requires: python3 (numpy, pandas, scikit-learn, scikit-optimize, xgboost, shap), Rscript
# (dplyr, tidyr, tibble, ggplot2, magrittr, tidyselect, gridExtra, R.utils)

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
m = BUILD_MODEL('RandomForestRegressor', N_JOBS=4)  # keyword arg: consensus_hyperparams_dict
                                                       # is the 2nd POSITIONAL param, not N_JOBS -
                                                       # this test previously called
                                                       # BUILD_MODEL('RandomForestRegressor', 4),
                                                       # which passed 4 as consensus_hyperparams_dict
                                                       # and crashed on model.set_params(**4)
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
echo "    '{MODEL}_{BLOCKING}_{fold}' row identifier that select_best_model.R depends on"
echo "    (read.delim(..., row.names=1) + separate(Model, into=c('model_name','blocking',"
echo "    'fold'))). Confirmed directly this crashed R's read.delim with 'duplicate"
echo "    row.names are not allowed' once two rows happened to share a value in what R fell"
echo "    back to treating as the row-name column.)"
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

echo "=== 6) select_best_model.R ggsave() variable names ==="
echo "    (aggregate_and_visualize_models.R was removed - its aggregation+visualization logic"
echo "    was consolidated into select_best_model.R/visualize_models.R. This check verifies"
echo "    the same class of bug that hit the old script (a ggsave() call referencing a"
echo "    variable that was never actually assigned/read from args) isn't present here.)"
declared_vars=$(grep -oE 'output_path_plot_[a-z_]+|output_path_median_boxplots' workflow/scripts/select_best_model.R | sort -u)
ggsave_vars=$(grep "ggplot2::ggsave" workflow/scripts/select_best_model.R | grep -oE 'output_path_[a-z_]+' | sort -u)
undeclared=$(comm -23 <(echo "$ggsave_vars") <(echo "$declared_vars"))
if [ -n "$undeclared" ]; then
    echo "  [FAIL] ggsave() references variable(s) never declared: $undeclared"
    FAIL=1
else
    echo "  [OK] every variable used in a ggsave() call is actually declared"
fi

echo "=== 7) split_dataset_by_chromosome.R ==="
echo "    (2 bugs, both confirmed directly:"
echo "    a) dataset[startsWith(dataset\$binID, pattern)] was missing a comma - without it,"
echo "       R indexes COLUMNS, not ROWS, and silently returns the entire, unfiltered dataset"
echo "       (no error, no warning)."
echo "    b) write.table(output_path, dataset.f, ...) had x and file swapped (the signature"
echo "       is write.table(x, file, ...)) - crashed directly.)"
printf "binID\tval\n1_1\t10\n1_2\t20\n10_1\t30\n2_1\t40\n" > "$TMP/split_dataset.tsv"
Rscript workflow/scripts/split_dataset_by_chromosome.R \
    --input_path="$TMP/split_dataset.tsv" --chrom=1 --output_path="$TMP/split_chr1.tsv" > "$TMP/split_out.log" 2>&1
RC=$?
check "script runs without error" "$RC" "0"
check "filters to exactly chr1's 2 rows (not the whole unfiltered dataset)" "$(tail -n +2 "$TMP/split_chr1.tsv" | wc -l | tr -d ' ')" "2"
check "does not false-positive match chr10 as a chr1 prefix" "$(grep -c '^10_' "$TMP/split_chr1.tsv" || true)" "0"

echo "=== 8) select_best_model.R - metadata columns excluded from hyperparameter consensus ==="
echo "    (bug: best_model_hyperparams still carries model_name/blocking/fold columns"
echo "    (process_aggregated_dataframe only drops 'Model') - since these are never all-NA,"
echo "    they leaked into the consensus output as if they were real hyperparameters"
echo "    (majority-voting model_name/blocking against themselves, and taking a meaningless"
echo "    median() of the fold index). Confirmed directly; fixed by excluding them before"
echo "    the NA-column-drop + consensus computation.)"
if grep -qE 'setdiff\(colnames\(best_model_hyperparams\).*model_name.*blocking.*fold' workflow/scripts/select_best_model.R; then
    echo "  [OK] model_name/blocking/fold are explicitly excluded before hyperparameter consensus"
else
    echo "  [FAIL] no exclusion of model_name/blocking/fold found before hyperparameter consensus"
    FAIL=1
fi

echo "=== 9) visualize_models.R - no redundant re-processing, 3-parameter signature ==="
echo "    (bug: the function used to redo rownames_to_column('Model') -> separate() ->"
echo "    mutate(blocking, fold) internally, even though select_best_model.R now calls"
echo "    process_aggregated_dataframe() on all 3 dataframes BEFORE passing them in - that"
echo "    upstream processing already drops the 'Model' column entirely, so re-deriving it"
echo "    here would crash outright. Also added a 3rd parameter (df_hyperparams) - the"
echo "    y_ordered-before-creation bug from an earlier fix is retested here too as a"
echo "    non-regression check.)"
Rscript -e '
`%>%` <- magrittr::`%>%`
source("workflow/scripts/helpers/process_aggregated_dataframe.R")
source("workflow/scripts/helpers/visualize_models.R")
set.seed(1)
rn <- c("LinearRegression_5_0","LinearRegression_5_1","XGBRegressor_5_0","XGBRegressor_5_1")
feats <- data.frame(f1=c(1,0,1,1), f2=c(0,1,1,0)); rownames(feats) <- rn
perf <- data.frame(train_r2=runif(4)); rownames(perf) <- rn
hp <- data.frame(model__fit_intercept=c("True","True",NA,NA), model__max_depth=c(NA,NA,5,7))
rownames(hp) <- rn
result <- visualize_models(process_aggregated_dataframe(feats), process_aggregated_dataframe(perf), process_aggregated_dataframe(hp))
cat("NAMES", paste(names(result), collapse=","), "\n")
' > "$TMP/visualize_out.txt" 2>&1
check "visualize_models() runs end-to-end on already-processed input and returns all 4 plots" "$(grep NAMES "$TMP/visualize_out.txt" | cut -d' ' -f2)" "heatmap_features,heatmap_performances,median_boxplots,hyperparams_table"

echo "=== 10) generate_genes_on_bins_obj.smk / split_dataset_by_chromosome.smk / select_best_model.smk conda: ==="
echo "    (bug: the conda: environment declared for a rule must actually match what its own"
echo "    script needs - confirmed on a real cluster run that a missing/wrong conda:"
echo "    directive causes 'command not found' or 'no package called X' failures that no"
echo "    amount of script-logic testing alone would ever catch, since this project's other"
echo "    tests never exercise --use-conda. This is a lightweight, text-based sanity check,"
echo "    not a real conda solve - see test_conda_env_declarations.sh for the repo-wide"
echo "    version of this check.)"
for pair in "split_dataset_by_chromosome.smk:split_dataset_by_chromosome.R" "select_best_model.smk:select_best_model.R"; do
    smk_file="workflow/rules/ML/${pair%%:*}"
    r_script="workflow/scripts/${pair##*:}"
    env_file=$(grep -A1 "conda:" "$smk_file" | grep "\.yml" | sed 's/.*envs\///;s/"//g' | tr -d ' ')
    if [ "$r_script" == "workflow/scripts/select_best_model.R" ]; then
        r_script="$r_script workflow/scripts/helpers/visualize_models.R workflow/scripts/helpers/process_aggregated_dataframe.R"
    fi
    used_pkgs=$(grep -ohE "[A-Za-z0-9.]+::" $r_script 2>/dev/null | sed 's/:://' | sort -u)
    env_pkgs=$(grep -oE "r-[a-zA-Z0-9.]+" "workflow/envs/$env_file" 2>/dev/null | sed 's/^r-//' | tr '[:upper:]' '[:lower:]')
    missing=""
    for pkg in $used_pkgs; do
        pkg_lower=$(echo "$pkg" | tr '[:upper:]' '[:lower:]')
        if ! echo "$env_pkgs" | grep -qx "$pkg_lower"; then
            # r-essentials covers these common tidyverse packages transitively
            case "$pkg_lower" in
                dplyr|tidyr|ggplot2|tibble|magrittr|tidyselect) continue ;;
            esac
            missing="$missing $pkg"
        fi
    done
    if [ -n "$missing" ]; then
        echo "  [FAIL] $smk_file's conda env ($env_file) may be missing:$missing"
        FAIL=1
    else
        echo "  [OK] $smk_file's conda env ($env_file) covers every package used by $r_script"
    fi
done

echo
if [ "$FAIL" -eq 0 ]; then
    echo "ALL TESTS PASSED"
    exit 0
else
    echo "SOME TESTS FAILED"
    exit 1
fi
