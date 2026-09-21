import numpy as np
import pandas as pd
import math as mh
import scipy as sp
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

from workflow.scripts.helpers.BUILD_MODEL import *
from workflow.scripts.helpers.GET_HYPERPARAM_SPACE import *
from workflow.scripts.helpers.SPLIT_DATASET import *
from workflow.scripts.helpers.CHECK_CHR_SORTING import *
from workflow.scripts.helpers.BLOCK_DATASET import *
from workflow.scripts.helpers.GENERATE_KFOLD_IDXS import *
from workflow.scripts.helpers.SUBSET_DATA_BY_IDXS import *
from workflow.scripts.helpers.OPTIMIZE_HYPERPARAMS import *
from workflow.scripts.helpers.EVALUATE_PERFORMANCES import *
from workflow.scripts.helpers.AGGREGATE_SELECTED_FEATURES import *

warnings.filterwarnings("ignore", message="Skipping features without any observed values")

parser = argparse.ArgumentParser()

parser.add_argument('--input_path_dataset', dest='input_path_dataset', required=True, type=str, help='')
parser.add_argument('--input_path_is_autocorrelated', dest='input_path_is_autocorrelated', required=True, type=str, help='')
parser.add_argument('--input_path_is_multicollinear', dest='input_path_is_multicollinear', required=True, type=str, help='')
parser.add_argument('--input_path_models', dest='input_path_models', required=True, type=str, help='')
parser.add_argument('--input_path_blockings', dest='input_path_blockings', required=True, type=str, help='')
parser.add_argument('--target_col', dest = 'target_col', required = True, type = str, help = '')
parser.add_argument('--selected_model', dest = 'selected_model', required = True, type = str, help = '')
parser.add_argument('--selected_blocking', dest = 'selected_blocking', required = True, type = str, help = '')
parser.add_argument('--outer_cv', dest = 'outer_cv', required = True, type = int, help = '')
parser.add_argument('--inner_cv', dest = 'inner_cv', required = True, type = int, help = '')
parser.add_argument('--test_size', dest = 'test_size', required = True, type = float, help = '')
parser.add_argument('--optimizer', dest = 'optimizer', required = True, type = str, help = '')
parser.add_argument('--n_iter', dest = 'n_iter', required = True, type = int, help = '')
parser.add_argument('--n_jobs', dest = 'n_jobs', required = True, type = int, help = '')
parser.add_argument('--model_n_jobs', dest = 'model_n_jobs', required = True, type = int, help = '')
parser.add_argument('--output_dir_models', dest = 'output_dir_models', required = True, type = str, help = '')
parser.add_argument('--output_path_warning', dest = 'output_path_warning', required = True, type = str, help = '')
parser.add_argument('--output_path_hyperparams', dest = 'output_path_hyperparams', required = True, type = str, help = '')
parser.add_argument('--output_path_selected_features', dest = 'output_path_selected_features', required = True, type = str, help = '')
parser.add_argument('--output_path_performances', dest = 'output_path_performances', required = True, type = str, help = '')

args = parser.parse_args()

DATASET = pd.read_csv(args.input_path_dataset, header = 0, sep = "\t")

TARGET_COL = args.target_col
TEST_SIZE = args.test_size
OUTER_CV = args.outer_cv
INNER_CV = args.inner_cv
OPTIMIZER = args.optimizer
N_ITER = args.n_iter
N_JOBS = args.n_jobs
MODEL_N_JOBS = args.model_n_jobs
SEED = 48047510 #!# TODO: list of seeds

METRICS_LIST = ["r2", "adj_r2", "explained_variance", "neg_mean_squared_error",
                "neg_root_mean_squared_error", "neg_mean_absolute_error", "neg_median_absolute_error",
                "neg_max_error", "d2_absolute_error_score"]

with open(args.input_path_is_autocorrelated) as f:
    line = f.readline().rstrip()
    IS_AUTOCORRELATED = (line == "TRUE")

with open(args.input_path_is_multicollinear) as f:
    line = f.readline().rstrip()
    IS_MULTICOLLINEAR = (line == "TRUE")

models_df = pd.read_csv(args.input_path_models, header = 0, sep = "\t")
MODELS = models_df.model.tolist()

if args.selected_model != "all":
    MODELS = [MODELS[MODELS.index(args.selected_model)]]

with open(args.input_path_blockings) as f:
    BLOCKINGS = [line.rstrip() for line in f.readlines()]
    if args.selected_blocking != "all":
        BLOCKINGS = [BLOCKINGS[BLOCKINGS.index(args.selected_blocking)]]

num_rows_dataset = DATASET.shape[0]

best_hyperparams_dict = dict()
selected_features_dict = dict()
performances_dict = dict()

for MODEL in MODELS:

    HANDLES_NAs = bool(models_df[models_df.model == MODEL].handles_NAs.values)
    NEEDS_SCALING = bool(models_df[models_df.model == MODEL].needs_scaling.values)
    IS_TREE_BASED = bool(models_df[models_df.model == MODEL].is_tree_based.values)
    IS_LINEAR = bool(models_df[models_df.model == MODEL].is_linear.values)
    REGRESSOR = BUILD_MODEL(MODEL, N_JOBS=MODEL_N_JOBS)
    HYPERPARAM_SPACE = GET_HYPERPARAM_SPACE(MODEL) # add hyperparam spaces
    bool_optimize_hyperparams = True

    for BLOCKING in BLOCKINGS:

        X_train, X_test, \
        y_train, y_test, \
        binIDs_train, binIDs_test = SPLIT_DATASET(
            DATASET,
            TARGET_COL = TARGET_COL,
            COLS_TO_EXCLUDE = ["FDR"],
            BLOCKING = BLOCKING,
            SEED = SEED,
            TEST_SIZE = TEST_SIZE,
            CHECK_SORTING_FUN = CHECK_CHR_SORTING,
            BLOCK_DATASET_FUN = BLOCK_DATASET
        )

        num_rows_test_set = X_test.shape[0]
        empirical_test_size = (num_rows_test_set / num_rows_dataset)

        outer_idxs = GENERATE_KFOLD_IDXS(BLOCKING, OUTER_CV, X_train, binIDs_train, CHECK_CHR_SORTING)

        for outer_fold, outer_idxs_tuple in outer_idxs.items():

            row_name = f"{MODEL}_{BLOCKING}_{outer_fold}"
            print(row_name)

            X_train_train, X_train_val, \
            y_train_train, y_train_val, \
            binIDs_train_train, binIDs_train_val = SUBSET_DATA_BY_IDXS(X_train, y_train, binIDs_train, outer_idxs_tuple)

            if (not NEEDS_SCALING) and (HANDLES_NAs):

                pipe = sk.pipeline.Pipeline([
                    ("feature_selection", sk.feature_selection.SequentialFeatureSelector(estimator = sk.base.clone(REGRESSOR), n_features_to_select = "auto", direction = 'forward', cv = INNER_CV, n_jobs = N_JOBS)),
                    ("model", sk.base.clone(REGRESSOR))
                ])

            else:

                pipe = sk.pipeline.Pipeline([
                    ("imputation", sk.impute.SimpleImputer()),
                    ("scaling", sk.preprocessing.StandardScaler()),
                    ("feature_selection", sk.feature_selection.SequentialFeatureSelector(estimator = sk.base.clone(REGRESSOR), n_features_to_select = "auto", direction = 'forward', cv = INNER_CV, n_jobs = N_JOBS)),
                    ("model", sk.base.clone(REGRESSOR))
                ])

            best_model, best_hyperparams = OPTIMIZE_HYPERPARAMS(pipe, HYPERPARAM_SPACE, OPTIMIZER, X_train_train, y_train_train, SEED, N_JOBS, N_ITER, INNER_CV)

            performances = EVALUATE_PERFORMANCES(best_model, X_train_val, y_train_val, 1, METRICS_LIST)

            best_hyperparams_dict[row_name] = best_hyperparams

            imputed_columns = X_train_train.columns
            if "imputation" in best_model.named_steps:
                imputed_columns = best_model.named_steps["imputation"].get_feature_names_out(X_train_train.columns)
            selected_features_dict[row_name] = pd.Index(imputed_columns)[best_model.named_steps["feature_selection"].get_support()]

            performances_dict[row_name] = performances

            jbl.dump(best_model, (args.output_dir_models + "/" + row_name + ".pkl"))

            with open(args.output_path_warning, "a") as f:
                f.write(f"[{outer_fold}] {args.input_path_dataset} with BLOCKING = {str(BLOCKING)} has empirical TEST_SIZE = {str(empirical_test_size)}" + "\n")
                f.write(f"[{outer_fold}] For {MODEL} bool_optimize_hyperparams = {str(bool_optimize_hyperparams)}" + "\n")

                for k,v in performances.items():
                    f.write(f"[{outer_fold}] {MODEL} scored {v} on {k} with BLOCKING = {str(BLOCKING)}" + "\n")

            if BLOCKING == "noblocking" and IS_AUTOCORRELATED:
                with open(args.output_path_warning, "a") as f:
                    f.write(f"[{outer_fold}] {MODEL} was trained on {args.input_path_dataset} with BLOCKING = {str(BLOCKING)} while IS_AUTOCORRELATED = {str(IS_AUTOCORRELATED)}" + "\n")


aggregated_best_hyperparams = pd.DataFrame.from_dict(best_hyperparams_dict, orient="index")
aggregated_selected_features = AGGREGATE_SELECTED_FEATURES(X_train.columns.tolist(), selected_features_dict)
aggregated_performances = pd.DataFrame.from_dict(performances_dict, orient="index")

aggregated_best_hyperparams.to_csv(args.output_path_hyperparams, sep="\t", header=True, index=True)
aggregated_selected_features.to_csv(args.output_path_selected_features, sep="\t", header=True, index=True)
aggregated_performances.to_csv(args.output_path_performances, sep="\t", header=True, index=True)