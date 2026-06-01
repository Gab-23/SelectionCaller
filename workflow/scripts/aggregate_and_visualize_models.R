source("workflow/scripts/helpers/visualize_models.R")

options(scipen=999)        

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

input_dir_performance                <-  args$input_dir_performance
input_dir_selfeats                   <-  args$input_dir_selfeats
input_dir_hyperparam                 <-  args$input_dir_hyperparam
output_path_performance              <-  args$output_path_performance
output_path_selfeats                 <-  args$output_path_selfeats
output_path_median_boxplots          <-  args$output_path_median_boxplots

`%>%` <- magrittr::`%>%`

performances_paths <- list.files(input_dir_performance, pattern = "perf_file", full.names = T, recursive = T)
selfeats_paths <- list.files(input_dir_selfeats, pattern = "selfeats_file", full.names = T, recursive = T)

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

performances_aggregated <- do.call(rbind, performances_list)
selfeats_aggregated <- do.call(rbind, selfeats_list)

output_visualization <- visualize_models(selfeats_aggregated, performances_aggregated)

ggplot2::ggsave(output_path_heatmap_sel_feats, output_visualization$heatmap_features, width = 20, height = 10)
ggplot2::ggsave(output_path_heatmap_performances, output_visualization$heatmap_performances, width = 15, height = 10)
ggplot2::ggsave(output_path_median_boxplots, output_visualization$median_boxplots, width = 15, height = 10)