source("workflow/scripts/helpers/visualize_models.R")
source("workflow/scripts/helpers/process_aggregated_dataframe.R")

options(scipen=999)                                                               # unable scientific notation  

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_dir_performances               <- args$input_dir_performances
input_dir_hyperparams                <- args$input_dir_hyperparams
input_dir_selfeats                   <- args$input_dir_selfeats

output_path_best_model_name          <- args$output_path_best_model_name
output_path_consensus_hyperparams    <- args$output_path_consensus_hyperparams
output_path_consensus_features       <- args$output_path_consensus_features
output_path_plot_selfeats            <- args$output_path_plot_selfeats
output_path_plot_performance         <- args$output_path_plot_performance
output_path_median_boxplots          <- args$output_path_median_boxplots
output_path_plot_hyperparams         <- args$output_path_plot_hyperparams

eval_criterion <- "train_r2"

`%>%` <- magrittr::`%>%`

performances_paths <- list.files(input_dir_performances, pattern = "perf_file", full.names = T, recursive = T)
selfeats_paths <- list.files(input_dir_selfeats, pattern = "selfeats_file", full.names = T, recursive = T)
hyperparams_path <- list.files(input_dir_hyperparams, pattern = "hyperparams", full.names = T, recursive = T)

performances_list <- lapply(X = performances_paths, 
                            FUN = function(path){
                                df <- read.delim(path, sep = "\t", row.names = 1, check.names = F)
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

performances_aggregated <- do.call(rbind, performances_list)
selfeats_aggregated <- do.call(rbind, selfeats_list)

hyperparams_aggregated <- hyperparams_list |>
  lapply(tibble::rownames_to_column, var = "Model") |>
  dplyr::bind_rows() |>
  tibble::column_to_rownames("Model")


performances_aggregated <- process_aggregated_dataframe(performances_aggregated)
selfeats_aggregated <- process_aggregated_dataframe(selfeats_aggregated)
hyperparams_aggregated <- process_aggregated_dataframe(hyperparams_aggregated)

output_visualization <- visualize_models(selfeats_aggregated, performances_aggregated, hyperparams_aggregated)

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
writeLines(con = output_path_consensus_features, text = most_common_features)
write.table(file = output_path_consensus_hyperparams, x = consensus_hyperparams, sep = "\t", row.names = F, col.names = F, quote = F)

ggplot2::ggsave(output_path_plot_selfeats, output_visualization$heatmap_features, width = 20, height = 10)
ggplot2::ggsave(output_path_plot_performance, output_visualization$heatmap_performances, width = 15, height = 10)
ggplot2::ggsave(output_path_median_boxplots, output_visualization$median_boxplots, width = 15, height = 10)
ggplot2::ggsave(output_path_plot_hyperparams, output_visualization$hyperparams_table, width = 15, height = 10)