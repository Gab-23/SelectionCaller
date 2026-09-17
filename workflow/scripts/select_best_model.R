source("workflow/scripts/helpers/visualize_models.R")
source("workflow/scripts/helpers/process_aggregated_dataframe.R")

options(scipen=999)                                                               # unable scientific notation  

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

results_dir                          <- args$results_dir
input_dir_hyperparams                <- args$input_dir_hyperparams
input_dir_selfeats                   <- args$input_dir_selfeats

output_path_best_model_name          <- args$output_path_best_model_name
output_path_best_blocking            <- args$output_path_best_blocking
output_path_consensus_hyperparams    <- args$output_path_consensus_hyperparams
output_path_consensus_features       <- args$output_path_consensus_features
output_path_plot_selfeats            <- args$output_path_plot_selfeats
output_path_plot_performance         <- args$output_path_plot_performance
output_path_median_boxplots          <- args$output_path_median_boxplots
output_path_plot_hyperparams         <- args$output_path_plot_hyperparams

eval_criterion <- "train_r2"

`%>%` <- magrittr::`%>%`

# NOTE (fix, pooled model selection): performances are now read from EVERY chromosome's own
# per_chrom/{chrom}/performances/ directory (results_dir/per_chrom/*/performances/), not just
# this job's own chromosome - the (model, blocking) choice is pooled genome-wide, so every
# chromosome's job reads the same full set and independently arrives at the same winner (see
# conversation: avoids confounding "difference due to the chromosome" with "difference due to
# model architecture" when later comparing feature importance across chromosomes).
# input_dir_hyperparams/input_dir_selfeats stay chromosome-specific: once the winning
# (model, blocking) is known, hyperparameter/feature consensus is computed from ONLY this
# chromosome's own folds - that is where the chromosome-specific signal needs to live.
performances_paths <- list.files(file.path(results_dir, "per_chrom"), pattern = "perf_file", full.names = T, recursive = T)
selfeats_paths <- list.files(input_dir_selfeats, pattern = "selfeats_file", full.names = T, recursive = T)
hyperparams_path <- list.files(input_dir_hyperparams, pattern = "hyperparams", full.names = T, recursive = T)

performances_list <- lapply(X = performances_paths, 
                            FUN = function(path){
                                df <- read.delim(path, sep = "\t", row.names = 1, check.names = F)
                                # NOTE (fix): rownames alone ("{model}_{blocking}_{fold}") are NOT
                                # unique once pooled across chromosomes - every chromosome
                                # independently produces the same row names for the same
                                # (model, blocking, fold), and do.call(rbind, ...) silently mangles
                                # duplicate ROWNAMES (confirmed directly: "XGBRegressor_5_0" from
                                # two different chromosomes becomes "XGBRegressor_5_0" and
                                # "XGBRegressor_5_01" - corrupting the fold number once
                                # re-parsed). A regular data-frame COLUMN, unlike rownames, has no
                                # such uniqueness constraint - duplicate values across rows are
                                # completely fine (confirmed directly) - so converting rownames to
                                # an explicit "Model" column before combining avoids the collision
                                # entirely, with no need to alter the value itself. An explicit
                                # "chrom" column (parsed from the file path:
                                # .../per_chrom/{chrom}/performances/perf_file_*.tsv) is also kept
                                # for provenance/traceability.
                                chrom <- basename(dirname(dirname(path)))
                                df <- tibble::rownames_to_column(df, "Model")
                                df$chrom <- chrom
                                return(df)
                    })

selfeats_list <- lapply(X = selfeats_paths, 
                        FUN = function(path){
                            df <- read.delim(path, sep = "\t", row.names = 1, check.names = F)
                            return(df)
                    })

hyperparams_list <- lapply(X = hyperparams_path, 
                        FUN = function(path){
                            df <- read.delim(path, sep = "\t", row.names = 1, check.names = F)
                            return(df)
                    })

performances_aggregated <- dplyr::bind_rows(performances_list)
selfeats_aggregated <- do.call(rbind, selfeats_list)

hyperparams_aggregated <- hyperparams_list |>
  lapply(tibble::rownames_to_column, var = "Model") |>
  dplyr::bind_rows() |>
  tibble::column_to_rownames("Model")


performances_aggregated <- process_aggregated_dataframe(performances_aggregated)
selfeats_aggregated <- process_aggregated_dataframe(selfeats_aggregated)
hyperparams_aggregated <- process_aggregated_dataframe(hyperparams_aggregated)

output_visualization <- visualize_models(selfeats_aggregated, dplyr::select(performances_aggregated, -chrom), hyperparams_aggregated)

mean_performance <- performances_aggregated %>% 
  dplyr::group_by(model_name, blocking) %>% 
  dplyr::summarise(
    !!paste0("avg_", eval_criterion) := mean(.data[[eval_criterion]], na.rm = TRUE),
    .groups = "drop"
  )

best_model_idx <- which.max(mean_performance$avg_train_r2)
best_model <- mean_performance[best_model_idx,]

best_model_name <- mean_performance$model_name[best_model_idx] # get best model in general
best_model_blocking <- mean_performance$blocking[best_model_idx]

# NOTE: selfeats_aggregated/hyperparams_aggregated are chromosome-specific (see above), so
# this restricts to the genome-wide winning (model, blocking) WITHIN this chromosome's own
# folds only.
best_model_selfeats <- selfeats_aggregated[(selfeats_aggregated$model_name == best_model_name & selfeats_aggregated$blocking == best_model_blocking), , drop = FALSE]
best_model_hyperparams <- hyperparams_aggregated[(hyperparams_aggregated$model_name == best_model_name & hyperparams_aggregated$blocking == best_model_blocking), , drop = FALSE]

# get most common features
feature_cols <- setdiff(colnames(best_model_selfeats), c("model_name", "blocking", "fold"))
most_common_features <- names(colSums(best_model_selfeats[, feature_cols, drop = FALSE])[colSums(best_model_selfeats[, feature_cols, drop = FALSE]) > ceiling(nrow(best_model_selfeats) / 2)])

# get median / more frequent hyperparams
# NOTE (fix): excluding model_name/blocking/fold here too, same reasoning as feature_cols above -
# best_model_hyperparams still carries these 3 metadata columns (process_aggregated_dataframe only
# drops "Model"), and since they're never all-NA they'd otherwise leak into the consensus output
# as if they were real hyperparameters (majority-voting model_name/blocking against themselves,
# and taking a meaningless median() of the fold index).
hyperparam_cols <- setdiff(colnames(best_model_hyperparams), c("model_name", "blocking", "fold"))
best_model_hyperparams <- best_model_hyperparams[, hyperparam_cols, drop = FALSE]

to.keep <- names(colSums(is.na(best_model_hyperparams))[colSums(is.na(best_model_hyperparams)) < nrow(best_model_hyperparams)])
best_model_hyperparams.s <- best_model_hyperparams[, to.keep, drop = FALSE]

consensus_hyperparams <- sapply(best_model_hyperparams.s, function(col) {
  if (is.numeric(col)) {
    as.character(median(col, na.rm = TRUE))
  } else {
    names(which.max(table(col)))
  }
})

consensus_hyperparams <- data.frame(
  hyperparam = names(consensus_hyperparams),
  median_hyperparam = unname(consensus_hyperparams),
  row.names = NULL
)

writeLines(con = output_path_best_model_name, text = best_model_name)
writeLines(con = output_path_best_blocking, text = as.character(best_model_blocking))
writeLines(con = output_path_consensus_features, text = most_common_features, sep = "\n") # remember that this is chromosome-specific
write.table(file = output_path_consensus_hyperparams, x = consensus_hyperparams, sep = "\t", row.names = F, col.names = F, quote = F)

ggplot2::ggsave(output_path_plot_selfeats, output_visualization$heatmap_features, width = 20, height = 10)
ggplot2::ggsave(output_path_plot_performance, output_visualization$heatmap_performances, width = 15, height = 10)
ggplot2::ggsave(output_path_median_boxplots, output_visualization$median_boxplots, width = 15, height = 10)
ggplot2::ggsave(output_path_plot_hyperparams, output_visualization$hyperparams_table, width = 15, height = 10)