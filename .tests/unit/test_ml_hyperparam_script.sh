#!/usr/bin/env bash
# Unit tests for how train_cv_plus_models.py reads the consensus hyperparameters written by
# select_best_model.R (workflow/scripts/helpers/coerce_hyperparam_value.py).
# Usage: bash .tests/unit/test_ml_hyperparams_scripts.sh
# Requires: python3 (pandas, numpy, scikit-learn, xgboost)

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
py() { PYTHONPATH="$REPO_ROOT" python3 -c "$1" 2>&1 | tail -1; }

echo "=== 1) coerce_hyperparam_value: every kind of value keeps its type ==="
echo "    (bug: values already parsed by pandas went straight into int(): True -> 1 and every"
echo "    fractional float was truncated, e.g. learning_rate 0.1 -> 0)"
IMPORT="from workflow.scripts.helpers.coerce_hyperparam_value import coerce_hyperparam_value as c; import numpy as np"
check "'TRUE' / 'FALSE' -> bool"       "$(py "$IMPORT; print([c('TRUE'), c('FALSE')], [type(c('TRUE')).__name__])")"  "[True, False] ['bool']"
check "'250' -> int, '0.1' -> float"   "$(py "$IMPORT; print(repr(c('250')), repr(c('0.1')))")"                      "250 0.1"
check "'rbf' -> str"                   "$(py "$IMPORT; print(repr(c('rbf')))")"                                       "'rbf'"
check "bool True / numpy bool -> bool (was int 1)" "$(py "$IMPORT; print(repr(c(True)), repr(c(np.bool_(False))))")"  "True False"
check "float 0.1 -> 0.1 (was truncated to 0)"      "$(py "$IMPORT; print(repr(c(0.1)), repr(c(np.float64(0.85))))")" "0.1 0.85"
check "whole-number float 4.0 -> int 4"            "$(py "$IMPORT; print(repr(c(4.0)), repr(c(np.int64(7))))")"     "4 7"

echo "=== 2) read_consensus_hyperparams on the two files that triggered the bug ==="
echo "    (LinearRegression: only booleans -> pandas made it a bool column, fit_intercept=1;"
echo "    XGBRegressor: only numbers -> pandas made it a float column, fractions truncated to 0)"
printf 'model__fit_intercept\tTRUE\nmodel__positive\tFALSE\n' > "$TMP/lr.tsv"
printf 'model__learning_rate\t0.1\nmodel__max_depth\t7\nmodel__subsample\t0.85\nmodel__colsample_bytree\t0.7\nmodel__gamma\t0.25\nmodel__reg_alpha\t0.001\nmodel__n_estimators\t250\n' > "$TMP/xgb.tsv"
READ="from workflow.scripts.helpers.coerce_hyperparam_value import read_consensus_hyperparams as r"
check "LinearRegression file -> real booleans, 'model__' prefix removed" \
    "$(py "$READ; d = r('$TMP/lr.tsv'); print(d, [type(v).__name__ for v in d.values()])")" \
    "{'fit_intercept': True, 'positive': False} ['bool', 'bool']"
check "XGBRegressor all-numeric file -> fractions kept, integers stay int" \
    "$(py "$READ; print(r('$TMP/xgb.tsv'))")" \
    "{'learning_rate': 0.1, 'max_depth': 7, 'subsample': 0.85, 'colsample_bytree': 0.7, 'gamma': 0.25, 'reg_alpha': 0.001, 'n_estimators': 250}"

echo "=== 3) the consensus models actually fit and learn (BUILD_MODEL + set_params) ==="
FIT="
import numpy as np
from workflow.scripts.helpers.coerce_hyperparam_value import read_consensus_hyperparams as r
from workflow.scripts.helpers.BUILD_MODEL import BUILD_MODEL
rng = np.random.default_rng(0); X = rng.normal(size=(300, 4)); y = 2 * X[:, 0] + rng.normal(scale=0.3, size=300)
"
check "LinearRegression with consensus hyperparameters fits (no 'fit_intercept ... Got 1')" \
    "$(py "$FIT
m = BUILD_MODEL('LinearRegression', r('$TMP/lr.tsv')).fit(X, y); print('fitted')")" "fitted"
check "XGBRegressor with consensus hyperparameters learns (was a constant prediction, R2 = 0)" \
    "$(py "$FIT
p = BUILD_MODEL('XGBRegressor', r('$TMP/xgb.tsv')).fit(X, y).predict(X)
print('learns' if p.std() > 0.5 and 1 - ((y - p)**2).sum() / ((y - y.mean())**2).sum() > 0.8 else 'constant')")" "learns"

echo "=== 4) train_cv_plus_models.py uses the shared reader (no local re-implementation) ==="
check "train_cv_plus_models.py calls read_consensus_hyperparams" \
    "$(grep -c 'read_consensus_hyperparams(args.input_path_consensus_hyperparams)' workflow/scripts/train_cv_plus_models.py)" "1"
check "train_cv_plus_models.py compiles" "$(python3 -m py_compile workflow/scripts/train_cv_plus_models.py && echo ok)" "ok"

echo
if [ "$FAIL" -eq 0 ]; then
    echo "ALL TESTS PASSED"
    exit 0
else
    echo "SOME TESTS FAILED"
    exit 1
fi