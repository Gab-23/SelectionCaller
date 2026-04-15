source("workflow/scripts/helpers/normATAC.R")
source("workflow/scripts/helpers/svaBatchCor.R")
source("workflow/scripts/helpers/pca_plot.R")
source("workflow/scripts/helpers/pca_scree_plot.R")
source("workflow/scripts/helpers/plot_distances_PCA.R")
source("workflow/scripts/helpers/do_get_counts_analysis.R")
source("workflow/scripts/helpers/plot_CNVkit_lines.R")

library(edgeR)
library(sva)

options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

#optional_analysis		                      <- FALSE
#correct_for_CNAs		                      <- as.logical(args$correct_for_CNAs)

input_dir_coverage                            <- args$input_dir_coverage
input_dir_scaling_factors	                  <- args$input_dir_scaling_factors
input_path_metadata                           <- args$input_path_metadata
input_path_samples_to_exclude                 <- args$input_path_samples_to_exclude

minSamples_atacseq                            <- as.integer(args$minSamples_atacseq)                  # TODO: this can be automated in the snakefile, based on the metadata file
minCounts_atacseq                             <- as.integer(args$minCounts_atacseq)                   # TODO: this can be automated in the snakefile, based on the metadata file
minCpm_atacseq                                <- as.integer(args$minCpm_atacseq)                      # TODO: this can be automated in the snakefile, based on the metadata file
SVs				                              <- as.integer(args$SVs)

output_path_library_sizes                     <- args$output_path_library_sizes
output_path_library_sizes_woCNAs              <- args$output_path_library_sizes_woCNAs

output_path_counts		                      <- args$output_path_counts
output_path_counts_woCNAs		              <- args$output_path_counts_woCNAs

#output_path_counts_norm                      <- args$output_path_counts_norm
#output_path_counts_sva                       <- args$output_path_counts_sva			      # this is just for plotting
output_path_model_matrix	                  <- args$output_path_model_matrix
output_path_model_matrix_woCNAs	              <- args$output_path_model_matrix_woCNAs

output_path_PCA_condition                     <- args$output_path_PCA_condition
output_path_PCA_condition_woCNAs              <- args$output_path_PCA_condition_woCNAs

output_path_PCA_library_size                  <- args$output_path_PCA_library_size
output_path_PCA_library_size_woCNAs           <- args$output_path_PCA_library_size_woCNAs

output_path_PCA_sex                           <- args$output_path_PCA_sex
output_path_PCA_sex_woCNAs                    <- args$output_path_PCA_sex_woCNAs

output_path_PCA_donor                         <- args$output_path_PCA_donor
output_path_PCA_donor_woCNAs                  <- args$output_path_PCA_donor_woCNAs

output_path_PCA_study                         <- args$output_path_PCA_study
output_path_PCA_study_woCNAs                  <- args$output_path_PCA_study_woCNAs

output_path_PCA_subtype                       <- args$output_path_PCA_subtype
output_path_PCA_subtype_woCNAs                <- args$output_path_PCA_subtype_woCNAs

output_path_scree_plot		                  <- args$output_path_scree_plot
output_path_scree_plot_woCNAs		          <- args$output_path_scree_plot_woCNAs

output_path_PCA_sva_condition                 <- args$output_path_PCA_sva_condition
output_path_PCA_sva_condition_woCNAs          <- args$output_path_PCA_sva_condition_woCNAs

output_path_PCA_sva_library_size              <- args$output_path_PCA_sva_library_size
output_path_PCA_sva_library_size_woCNAs       <- args$output_path_PCA_sva_library_size_woCNAs

output_path_PCA_sva_sex                       <- args$output_path_PCA_sva_sex
output_path_PCA_sva_sex_woCNAs                <- args$output_path_PCA_sva_sex_woCNAs

output_path_PCA_sva_donor                     <- args$output_path_PCA_sva_donor
output_path_PCA_sva_donor_woCNAs              <- args$output_path_PCA_sva_donor_woCNAs

output_path_PCA_sva_study                     <- args$output_path_PCA_sva_study
output_path_PCA_sva_study_woCNAs              <- args$output_path_PCA_sva_study_woCNAs

output_path_PCA_sva_subtype                   <- args$output_path_PCA_sva_subtype 
output_path_PCA_sva_subtype_woCNAs            <- args$output_path_PCA_sva_subtype_woCNAs 

output_path_scree_plot_sva	                  <- args$output_path_scree_plot_sva
output_path_scree_plot_sva_woCNAs	          <- args$output_path_scree_plot_sva_woCNAs

output_path_distance_plot	                  <- args$output_path_distance_plot
output_path_distance_plot_woCNAs	          <- args$output_path_distance_plot_woCNAs

output_path_CNVkit_lineplot                   <- args$output_path_CNVkit_lineplot

files_coverage <- list.files(input_dir_coverage, recursive = T, full.names = T)
files_scaling_factors <- list.files(input_dir_scaling_factors, recursive = T, full.names = T, pattern = "scaling")

metadata <- read.table(input_path_metadata, header = T, sep = "\t")
samples_to_exclude <- readLines(input_path_samples_to_exclude)

sample_names_coverage <- unlist(lapply(X = files_coverage, FUN = function(x){strsplit(x, split = "/")[[1]][[5]]}))
sample_names_scaling_factors <- unlist(lapply(X = files_scaling_factors, FUN = function(x){strsplit(x, split = "/")[[1]][[5]]}))

cov_list <- lapply(X = seq(length(files_coverage)), 
                   FUN = function(idx){
                    file_coverage <- files_coverage[idx]
                    sample_name_coverage <- sample_names_coverage[idx]
                    tab_coverage <- read.table(file_coverage, header = F, sep = "\t")
                    colnames(tab_coverage) <- c("chr", "start", "end", "binID", sample_name_coverage)
                    to.return <- tab_coverage[,c("binID", sample_name_coverage)]
                    return(to.return)
                    })

scaling_factors_list <- lapply(X = seq(length(files_scaling_factors)),
                   FUN = function(idx){
                    file_scaling_factors <- files_scaling_factors[idx]
                    sample_name_scaling_factors <- sample_names_scaling_factors[idx]
                    tab_scaling_factors <- read.table(file_scaling_factors, header = F, sep = "\t")
                    colnames(tab_scaling_factors) <- c("chr", "start", "end", "binID", "segmentID", sample_name_scaling_factors)
                    to.return <- tab_scaling_factors[,c("binID", sample_name_scaling_factors)]
                    return(to.return)
                    })

counts <- Reduce(f = function(x,y){merge(x,y, by = "binID")}, x = cov_list)
counts <- counts[,!(colnames(counts) %in% samples_to_exclude)]

scaling_factors <- Reduce(f = function(x, y){merge(x,y, by = "binID")}, x = scaling_factors_list)
scaling_factors <- scaling_factors[,!(colnames(scaling_factors) %in% samples_to_exclude)]

CNVkit_lineplot <- plot_CNVkit_lines(scaling_factors)

## filter counts to account bins in scaling_factors

out.list <- do_get_counts_analysis(counts = counts, scaling_factors = scaling_factors, correct_for_CNAs = FALSE)
out.list.woCNAs <- do_get_counts_analysis(counts = counts, scaling_factors = scaling_factors, correct_for_CNAs = TRUE)

############ FINAL STEPS

#if (optional_analysis){
#samples_to_keep_optional <- c("SRR24401466", "SRR24401465", "SRR27726641", "SRR27726640",
#                               "SRR24401468", "SRR24401467", "SRR24401462", "SRR24401461",
#                               "SRR27726632", "SRR27726631", "SRR27726628", "SRR27726627",
#                               "SRR27726634", "SRR27726633", "SRR27726630", "SRR27726629")
#
#conditiontumor_new <- ifelse(rownames(mm) %in% samples_to_keep_optional, ifelse(mm[,"conditiontumor"] == 1, 1, -1), 0)
#mm[,"conditiontumor"] <- conditiontumor_new
#
#}

############ SAVE OUTPUTS

write.table(x = out.list$count.f, file = output_path_counts, quote = F, sep = "\t", row.names = F, col.names = T)
#write.table(x = count.f.norm, file = output_path_counts_norm, quote = F, sep = "\t", row.names = F, col.names = T)
#write.table(x = count.f.norm.sva, file = output_path_counts_sva, quote = F, sep = "\t", row.names = F, col.names = T)
write.table(x = out.list$library_sizes_df, file = output_path_library_sizes, quote = F, sep = "\t", row.names = F, col.names = T)

save(out.list$mm, file = output_path_model_matrix)

ggplot2::ggsave(output_path_PCA_condition, plot = out.list$pca_plot_condition, width = 6, height = 5)
ggplot2::ggsave(output_path_PCA_library_size, plot = out.list$pca_plot_library_size, width = 6, height = 5)
ggplot2::ggsave(output_path_PCA_sex, plot = out.list$pca_plot_sex, width = 6, height = 5)
ggplot2::ggsave(output_path_PCA_donor, plot = out.list$pca_plot_donor, width = 6, height = 5)
ggplot2::ggsave(output_path_PCA_study, plot = out.list$pca_plot_study, width = 6, height = 5)
ggplot2::ggsave(output_path_PCA_subtype, plot = out.list$pca_plot_subtype, width = 6, height = 5)
ggplot2::ggsave(output_path_scree_plot, plot = out.list$scree_plot, width = 6, height = 5)

ggplot2::ggsave(output_path_PCA_sva_condition, plot = out.list$pca_sva_plot_condition, width = 6, height = 5)
ggplot2::ggsave(output_path_PCA_sva_library_size, plot = out.list$pca_sva_plot_library_size, width = 6, height = 5)
ggplot2::ggsave(output_path_PCA_sva_sex, plot = out.list$pca_sva_plot_sex, width = 6, height = 5)
ggplot2::ggsave(output_path_PCA_sva_donor, plot = out.list$pca_sva_plot_donor, width = 6, height = 5)
ggplot2::ggsave(output_path_PCA_sva_study, plot = out.list$pca_sva_plot_study, width = 6, height = 5)
ggplot2::ggsave(output_path_PCA_sva_subtype, plot = out.list$pca_sva_plot_subtype, width = 6, height = 5)
ggplot2::ggsave(output_path_scree_plot_sva, plot = out.list$scree_plot_sva, width = 6, height = 5)

ggplot2::ggsave(output_path_distance_plot, plot = out.list$pca_distance_plot, width = 6, height = 5)

write.table(x = out.list.woCNAs$count.f, file = output_path_counts, quote = F, sep = "\t", row.names = F, col.names = T)
#write.table(x = count.f.norm, file = output_path_counts_norm, quote = F, sep = "\t", row.names = F, col.names = T)
#write.table(x = count.f.norm.sva, file = output_path_counts_sva, quote = F, sep = "\t", row.names = F, col.names = T)
write.table(x = out.list.woCNAs$library_sizes_df, file = output_path_library_sizes, quote = F, sep = "\t", row.names = F, col.names = T)

save(out.list.woCNAs$mm, file = output_path_model_matrix)

ggplot2::ggsave(output_path_PCA_condition, plot = out.list.woCNAs$pca_plot_condition, width = 6, height = 5)
ggplot2::ggsave(output_path_PCA_library_size, plot = out.list.woCNAs$pca_plot_library_size, width = 6, height = 5)
ggplot2::ggsave(output_path_PCA_sex, plot = out.list.woCNAs$pca_plot_sex, width = 6, height = 5)
ggplot2::ggsave(output_path_PCA_donor, plot = out.list.woCNAs$pca_plot_donor, width = 6, height = 5)
ggplot2::ggsave(output_path_PCA_study, plot = out.list.woCNAs$pca_plot_study, width = 6, height = 5)
ggplot2::ggsave(output_path_PCA_subtype, plot = out.list.woCNAs$pca_plot_subtype, width = 6, height = 5)
ggplot2::ggsave(output_path_scree_plot, plot = out.list.woCNAs$scree_plot, width = 6, height = 5)

ggplot2::ggsave(output_path_PCA_sva_condition, plot = out.list.woCNAs$pca_sva_plot_condition, width = 6, height = 5)
ggplot2::ggsave(output_path_PCA_sva_library_size, plot = out.list.woCNAs$pca_sva_plot_library_size, width = 6, height = 5)
ggplot2::ggsave(output_path_PCA_sva_sex, plot = out.list.woCNAs$pca_sva_plot_sex, width = 6, height = 5)
ggplot2::ggsave(output_path_PCA_sva_donor, plot = out.list.woCNAs$pca_sva_plot_donor, width = 6, height = 5)
ggplot2::ggsave(output_path_PCA_sva_study, plot = out.list.woCNAs$pca_sva_plot_study, width = 6, height = 5)
ggplot2::ggsave(output_path_PCA_sva_subtype, plot = out.list.woCNAs$pca_sva_plot_subtype, width = 6, height = 5)
ggplot2::ggsave(output_path_scree_plot_sva, plot = out.list.woCNAs$scree_plot_sva, width = 6, height = 5)

ggplot2::ggsave(output_path_distance_plot, plot = out.list.woCNAs$pca_distance_plot, width = 6, height = 5)

ggplot2::ggsave(output_path_CNVkit_lineplot, plot = CNVkit_lineplot, width = 14, height = 4)
