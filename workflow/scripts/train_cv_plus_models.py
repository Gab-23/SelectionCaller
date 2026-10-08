import numpy as np
import pandas as pd
import math as mh
import scipy as sp
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
import sklearn as sk
import skopt as sko
import xgboost as xgb
import shap as sh
import pickle as pkl
import joblib as jbl
import seaborn as sns
import matplotlib.pyplot as plt
import random
import warnings
import argparse
import os

from workflow.scripts.helpers.coerce_hyperparam_value import *
from workflow.scripts.helpers.BUILD_MODEL import *
from workflow.scripts.helpers.SPLIT_DATASET import *
from workflow.scripts.helpers.CHECK_CHR_SORTING import *
from workflow.scripts.helpers.GENERATE_KFOLD_IDXS import *
from workflow.scripts.helpers.SUBSET_DATA_BY_IDXS import *
from workflow.scripts.helpers.EVALUATE_PERFORMANCES import *
from workflow.scripts.helpers.GET_OUT_OF_FOLD_PREDICTIONS import *
from workflow.scripts.helpers.GET_EXPLAINER import *
from workflow.scripts.helpers.GET_OUT_OF_FOLD_SHAP_VALUES import *
from workflow.scripts.helpers.VISUALIZE_SHAP_VALUES import *
from workflow.scripts.helpers.COMPUTE_MAD_Z_PVALS import *
from workflow.scripts.helpers.CHROMOSOME_FOLDS import CHROMOSOME_BLOCKING, CHROMOSOME_KFOLD_IDXS

parser = argparse.ArgumentParser()

parser.add_argument('--input_path_dataset', dest = 'input_path_dataset', required = True, type = str, help = '')
parser.add_argument('--input_path_models', dest = 'input_path_models', required = True, type = str, help = '')
parser.add_argument('--input_path_best_model_name', dest = 'input_path_best_model_name', required = True, type = str, help = '')
parser.add_argument('--input_path_best_blocking', dest = 'input_path_best_blocking', required = True, type = str, help = '')
parser.add_argument('--input_path_consensus_hyperparams', dest = 'input_path_consensus_hyperparams', required = True, type = str, help = '')
parser.add_argument('--input_path_consensus_features', dest = 'input_path_consensus_features', required = True, type = str, help = '')
parser.add_argument('--target_col', dest = 'target_col', required = True, type = str, help = '')
parser.add_argument('--fdr_threshold', dest = 'fdr_threshold', required = True, type = float, help = '')
parser.add_argument('--outer_cv', dest = 'outer_cv', required = True, type = int, help = '')
parser.add_argument('--model_n_jobs', dest = 'model_n_jobs', required = True, type = int, help = '')
parser.add_argument('--bin_size_bp', dest = 'bin_size_bp', required = True, type = int, help = '')
parser.add_argument('--output_dir_models', dest = 'output_dir_models', required = True, type = str, help = '')
parser.add_argument('--output_path_shap_values', dest = 'output_path_shap_values', required = True, type = str, help = '')
parser.add_argument('--output_path_shap_plots', dest = 'output_path_shap_plots', required = True, type = str, help = '')
parser.add_argument('--output_path_base_value', dest = 'output_path_base_value', required = True, type = str, help = '')
parser.add_argument('--output_path_significance_calls', dest = 'output_path_significance_calls', required = True, type = str, help = '')
parser.add_argument('--shap_plot_rolling_window_bins', dest = 'shap_plot_rolling_window_bins', required = False, type = int, default = 50, help = '')
parser.add_argument('--output_path_warning', dest = 'output_path_warning', required = True, type = str, help = '')

args = parser.parse_args()

DATASET = pd.read_csv(args.input_path_dataset, header = 0, sep = "\t")

METRICS_LIST = ["r2", "adj_r2", "explained_variance", "neg_mean_squared_error",
                "neg_root_mean_squared_error", "neg_mean_absolute_error", "neg_median_absolute_error",
                "neg_max_error", "d2_absolute_error_score"]

models_df = pd.read_csv(args.input_path_models, header = 0, sep = "\t")

# read as text and typed by coerce_hyperparam_value (see read_consensus_hyperparams for why)
consensus_hyperparams_dict = read_consensus_hyperparams(args.input_path_consensus_hyperparams)

with open(args.input_path_consensus_features) as f:
    CONSENSUS_FEATURES = [line.rstrip() for line in f.readlines()]

non_consensus_features = [c for c in DATASET.columns if c not in CONSENSUS_FEATURES]

with open(args.input_path_best_model_name) as f:
    BEST_MODEL_NAME = f.readline().rstrip()

with open(args.input_path_best_blocking) as f:
    BEST_BLOCKING = f.readline().rstrip()

TARGET_COL = args.target_col
FDR_THRESHOLD = args.fdr_threshold
OUTER_CV = args.outer_cv
MODEL_N_JOBS = args.model_n_jobs

num_rows_dataset = DATASET.shape[0]

HANDLES_NAs = bool(models_df[models_df.model == BEST_MODEL_NAME].handles_NAs.values)
NEEDS_SCALING = bool(models_df[models_df.model == BEST_MODEL_NAME].needs_scaling.values)
IS_TREE_BASED = bool(models_df[models_df.model == BEST_MODEL_NAME].is_tree_based.values)
IS_LINEAR = bool(models_df[models_df.model == BEST_MODEL_NAME].is_linear.values)
REGRESSOR = BUILD_MODEL(BEST_MODEL_NAME, consensus_hyperparams_dict, MODEL_N_JOBS)

X, y, binIDs = SPLIT_DATASET(
                    DATASET,
                    TARGET_COL = TARGET_COL,
                    COLS_TO_EXCLUDE = non_consensus_features,
                    EXIT_EARLY = True
                )

# NOTE (LOCO): BEST_BLOCKING == "chromosome" -> folds made of whole chromosomes
if BEST_BLOCKING == CHROMOSOME_BLOCKING:
    outer_idxs = CHROMOSOME_KFOLD_IDXS(binIDs, OUTER_CV)
else:
    outer_idxs = GENERATE_KFOLD_IDXS(BEST_BLOCKING, OUTER_CV, X, binIDs, CHECK_CHR_SORTING)

os.makedirs(args.output_dir_models, exist_ok=True)

out_of_fold_predictions_list = []
out_of_fold_shap_values_list = []
base_values_list = []

for outer_fold, outer_idxs_tuple in outer_idxs.items():

    row_name = f"{BEST_MODEL_NAME}_{BEST_BLOCKING}_{outer_fold}"
    print(row_name)

    X_train, X_val, \
    y_train, y_val, \
    binIDs_train, binIDs_val = SUBSET_DATA_BY_IDXS(X, y, binIDs, outer_idxs_tuple)

    if (not NEEDS_SCALING) and (HANDLES_NAs):

        pipe = sk.pipeline.Pipeline([("model", sk.base.clone(REGRESSOR))])

    else:

        pipe = sk.pipeline.Pipeline([
            ("imputation", sk.impute.SimpleImputer()),
            ("scaling", sk.preprocessing.StandardScaler()),
            ("model", sk.base.clone(REGRESSOR))
        ])

    # NOTE (fix): fit / evaluation / dump were indented inside the `else:` above, so for models
    # that need no scaling and handle NAs (e.g. XGBRegressor) the pipeline was never fitted.
    # Now run for both branches.
    pipe.fit(X_train, y_train) # fit

    performances = EVALUATE_PERFORMANCES(pipe, X_val, y_val, 1, METRICS_LIST)
    with open(args.output_path_warning, "a") as f:
        for k, v in performances.items():
            f.write(f"[{outer_fold}] {BEST_MODEL_NAME} scored {v} on {k} with BLOCKING = {str(BEST_BLOCKING)}\n")

    jbl.dump(pipe, os.path.join(args.output_dir_models, f"{outer_fold}.pkl"))

    out_of_fold_predictions_elem = GET_OUT_OF_FOLD_PREDICTIONS(pipe, X_val, binIDs_val)
    out_of_fold_predictions_list.append(out_of_fold_predictions_elem)

    EXPLAINER = GET_EXPLAINER(pipe, IS_TREE_BASED, IS_LINEAR, X_train)
    out_of_fold_shap_values_elem = GET_OUT_OF_FOLD_SHAP_VALUES(EXPLAINER, X_val, binIDs_val)
    out_of_fold_shap_values_list.append(out_of_fold_shap_values_elem)

    base_values_list.append({"outer_fold": outer_fold, "base_value": EXPLAINER.expected_value})


out_of_fold_predictions = pd.concat(out_of_fold_predictions_list, axis=0)
out_of_fold_shap_values = pd.concat(out_of_fold_shap_values_list, axis=0)

assert out_of_fold_predictions.index.is_unique
assert set(out_of_fold_predictions.index) == set(DATASET["binID"])
assert out_of_fold_shap_values.index.is_unique
assert set(out_of_fold_shap_values.index) == set(DATASET["binID"])
assert len(out_of_fold_predictions) == num_rows_dataset
assert len(out_of_fold_shap_values) == num_rows_dataset

observation = DATASET.set_index("binID").loc[:, TARGET_COL]
out_of_fold_residuals = (out_of_fold_predictions["prediction"] - observation).reindex(out_of_fold_predictions.index)

p_values, z_robust = COMPUTE_MAD_Z_PVALS(out_of_fold_residuals)
_, q_values, _, _ = multipletests(p_values, method="fdr_bh")

significance_calls = pd.DataFrame({
    "observation": observation.reindex(out_of_fold_predictions.index),
    "prediction": out_of_fold_predictions["prediction"],
    "residual": out_of_fold_residuals,
    "z_robust": z_robust,
    "p_value": p_values,
    "q_value": q_values,
    "significant": q_values < FDR_THRESHOLD,
}, index=out_of_fold_predictions.index)

RAW_VALUES = DATASET.set_index("binID")

VISUALIZE_SHAP_VALUES(
    out_of_fold_shap_values,
    RAW_VALUES,
    args.bin_size_bp,
    args.output_path_shap_plots,
    ROLLING_WINDOW_BINS=args.shap_plot_rolling_window_bins
)

out_of_fold_shap_values.to_csv(args.output_path_shap_values, sep="\t", header=True, index=True)
pd.DataFrame(base_values_list).to_csv(args.output_path_base_value, sep="\t", header=True, index=False)
significance_calls.to_csv(args.output_path_significance_calls, sep="\t", header=True, index=True)