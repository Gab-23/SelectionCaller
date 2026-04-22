import numpy as np
import pandas as pd
import math as mh
import scipy as sp
import sklearn as sk
import skopt as sko
import xgboost as xgb
import shap as sh
import pickle as pkl
import seaborn as sns
import matplotlib.pyplot as plt
import pickleshare
import random
import warnings
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--input_path_dataset', dest='input_path_dataset', required=True, type=str, help='')
parser.add_argument('--input_path_is_autocorrelated', dest='input_path_is_autocorrelated', required=True, type=str, help='')
parser.add_argument('--input_path_is_multicollinear', dest='input_path_dataset', required=True, type=str, help='')
parser.add_argument('--input_path_models', dest='input_path_models', required=True, type=str, help='')
parser.add_argument('--input_path_blockings', dest='input_path_blockings', required=True, type=str, help='')
parser.add_argument('--input_selected_model', dest = 'input_selected_model', required = True, type = str, help = '')
parser.add_argument('--outer_cv', dest = 'outer_cv', required = True, type = int, help = '')
parser.add_argument('--inner_cv', dest = 'inner_cv', required = True, type = int, help = '')
parser.add_argument('--test_size', dest = 'test_size', required = True, type = float, help = '')
parser.add_argument('--optimizer', dest = 'optimizer', required = True, type = str, help = '')
parser.add_argument('--n_iter', dest = 'n_iter', required = True, type = int, help = '')
parser.add_argument('--n_jobs', dest = 'n_jobs', required = True, type = int, help = '')
parser.add_argument('--output_path_warning', dest = 'output_path_warning', required = True, type = str, help = '')
parser.add_argument('--output_path_hyperparams_first_round', dest = 'output_path_hyperparams_first_round', required = True, type = str, help = '')
parser.add_argument('--output_path_hyperparams_second_round', dest = 'output_path_hyperparams_second_round', required = True, type = str, help = '')
parser.add_argument('--output_path_selected_features', dest = 'output_path_selected_features', required = True, type = str, help = '')
parser.add_argument('--output_path_performances', dest = 'output_path_performances', required = True, type = str, help = '')

args = parser.parse_args()

DATASET = pd.read_csv(args.input_path_dataset, header = 0, sep = "\t") #!# we iterate over datasets
TEST_SIZE = args.test_size
OUTER_CV = args.outer_cv
INNER_CV = args.inner_cv
OPTIMIZER = args.optimizer
N_ITER = args.n_iter
N_JOBS = args.n_jobs
SEED = 48047510 #!# TODO: list of seeds

METRICS_LIST = ["r2", "explained_variance", "neg_mean_squared_error",
                "neg_root_mean_squared_error", "neg_mean_absolute_error", "neg_median_absolute_error",
                "neg_mean_absolute_percentage_error", "neg_mean_squared_log_error", "neg_root_mean_squared_log_error",
                "neg_max_error", "neg_mean_poisson_deviance", "neg_mean_gamma_deviance", "d2_absolute_error_score"]

with open(args.input_path_is_autocorrelated) as f:
    line = f.readline()
    IS_AUTOCORRELATED = np.where(line == "TRUE", True, False).tolist()

with open(args.input_path_is_multicollinear) as f:
    line = f.readline()
    IS_MULTICOLLINEAR = np.where(line == "TRUE", True, False).tolist()

with open(args.input_path_models) as f:
    MODELS = [line.rstrip() for line in f.readlines()]
    if args.input_path_models != "all":
        MODELS = [MODELS[MODELS.index(args.input_path_models)]]

with open(args.input_path_blockings) as f:
    BLOCKINGS = [line.rstrip() for line in f.readlines()]
    if args.input_path_blockings != "all":
        BLOCKINGS = [BLOCKINGS[BLOCKINGS.index(args.input_path_blockings)]]

# with open(args.input_path_seeds) as f:
#   SEEDS = [line.rstrip() for line in f.readlines()]

num_rows_dataset = DATASET.shape[0]

best_hyperparams_first_round_dict = dict()
best_hyperparams_second_round_dict = dict()
selected_features_dict = dict()
performances_dict = dict()

for MODEL in MODELS:
    REGRESSOR = BUILD_MODEL(MODEL)
    HYPERPARAM_SPACE = GET_HYPERPARAM_SPACE(MODEL)
    bool_optimize_hyperparams = True if len(HYPERPARAM_SPACE) > 0 else False
    for BLOCKING in BLOCKINGS:

        X_train, X_test, y_train, y_test, binIDs_train, binIDs_test = SPLIT_DATASET(DATASET, BLOCKING, SEED, TEST_SIZE)

        num_rows_test_set = X_test.shape[0]
        empirical_test_size = (num_rows_test_set / num_rows_dataset)
    
       outer_idxs = GENERATE_KFOLD_IDXS(BLOCKING, OUTER_CV, X_train)
       
       for outer_fold in int(outer_idxs.keys()):

    #       generate datasets
            inner_idxs = GENERATE_KFOLD_IDXS(BLOCKING, OUTER_CV, X_train_train)

            for inner_fold in int(inner_idxs.keys()):

    #           generate datasets

                row_name = f"{MODEL}_{BLOCKING}_{outer_fold}_{inner_fold}"

                if bool_optimize_hyperparams:
                    REGRESSOR_I, best_hyperparams_first_round = OPTIMIZE_HYPERPARAMS(REGRESSOR, HYPERPARAM_SPACE, OPTIMIZER, X_train, y_train, CV, SEED, N_JOBS, N_ITER)
                else:
                    best_hyperparams_first_round = {}
                    REGRESSOR_I = REGRESSOR

                selected_features = SELECT_FEATURES(REGRESSOR_I, X_train, y_train, args.cv) #!# TODO
                X_train_reduced = X_train.loc[:, X_train.columns in selected_features]
                X_test_reduced = X_test.loc[:, X_test.columns in selected_features]

                if bool_optimize_hyperparams:
                    REGRESSOR_II, best_hyperparams_second_round = OPTIMIZE_HYPERPARAMS(REGRESSOR_I, HYPERPARAM_SPACE, OPTIMIZER, X_train_reduced, y_train, CV, SEED, N_JOBS, N_ITER)
                else:
                    best_hyperparams_second_round = {}
                    REGRESSOR_II = REGRESSOR

                performances = EVALUATE_PERFORMANCES(REGRESSOR_II, [X_train_reduced, X_test_reduced], [y_train, y_test], METRICS_LIST, CV)

                best_hyperparams_first_round_dict[row_name] = best_hyperparams_first_round
                best_hyperparams_second_round_dict[row_name] = best_hyperparams_second_round
                selected_features_dict[row_name] = selected_features
                performances_dict[row_name] = performances

                with open(args.output_path_warning, "a") as f:
                    f.write(f"[{outer_idx}:{inner_idx}] {args.input_path_dataset} with BLOCKING = {str(BLOCKING)} has empirical TEST_SIZE = {str(empirical_test_size)}" + "\n")
                    f.write(f"[{outer_idx}:{inner_idx}] For {MODEL} bool_optimize_hyperparams = {str(bool_optimize_hyperparams)}" + "\n")

                    for k,v in performances.items():
                        f.write(f"[{outer_idx}:{inner_idx}] {MODEL} scored {v} on {k} (CV = {CV}) with BLOCKING = {str(BLOCKING)}" + "\n")

                if BLOCKING == "no_blocking" and IS_AUTOCORRELATED:
                    with open(args.output_path_warning, "a") as f:
                        f.write(f"[{outer_idx}:{inner_idx}] {MODEL} was trained on {args.input_path_dataset} with BLOCKING = {str(BLOCKING)} while IS_AUTOCORRELATED = {str(IS_AUTOCORRELATED)}" + "\n")
    
        
aggregated_best_hyperparams_first_round = pd.DataFrame.from_dict(best_hyperparams_first_round_dict, orient="index")
aggregated_best_hyperparams_second_round = pd.DataFrame.from_dict(best_hyperparams_second_round_dict, orient="index")
aggregated_selected_features = AGGREGATE_SELECTED_FEATURES(selected_features_dict)
aggregated_performances = pd.DataFrame.from_dict(performances_dict, orient="index")

aggregated_best_hyperparams_first_round.to_csv(args.output_path_hyperparams_first_round, sep = "\t", header = True, index = False)
aggregated_best_hyperparams_second_round.to_csv(args.output_path_hyperparams_second_round, sep = "\t", header = True, index = False)
aggregated_selected_features.to_csv(args.output_path_selected_features, sep = "\t", header = True, index = False)
aggregated_performances.to_csv(args.output_path_performances, sep = "\t", header = True, index = False)
