source("workflow/scripts/helpers/edg2ATAC.R")
source("workflow/scripts/helpers/plot_diagnostics.R")
source("workflow/scripts/helpers/plot_volcano.R")
source("workflow/scripts/helpers/plot_landscape_FC.R")
#source("workflow/scripts/helpers/compute_FC_manually.R")
source("workflow/scripts/helpers/process_topTags.R")
source("workflow/scripts/helpers/plot_comparison_FC.R")
source("workflow/scripts/helpers/compare_acc_CNAs.R")
source("workflow/scripts/helpers/do_DA_analysis.R")

library(edgeR)
library(sva)

options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_metadata                            <- args$input_path_metadata
input_path_counts                              <- args$input_path_counts
input_path_counts_woCNAs                       <- args$input_path_counts_woCNAs
input_path_library_sizes                       <- args$input_path_library_sizes
input_path_library_sizes_woCNAs                <- args$input_path_library_sizes_woCNAs
input_path_model_matrix                        <- args$input_path_model_matrix
input_path_model_matrix_woCNAs                 <- args$input_path_model_matrix_woCNAs
input_path_ampl			                       <- args$input_path_ampl
input_path_del                                 <- args$input_path_del

output_path_DARs		                       <- args$output_path_DARs
output_path_DARs_woCNAs		                   <- args$output_path_DARs_woCNAs
output_path_diagnostics_plot	               <- args$output_path_diagnostics_plot
output_path_diagnostics_plot_woCNAs	           <- args$output_path_diagnostics_plot_woCNAs
output_path_volcano_plot                       <- args$output_path_volcano_plot
output_path_volcano_plot_woCNAs                <- args$output_path_volcano_plot_woCNAs
output_path_landscape_plot                     <- args$output_path_landscape_plot
output_path_landscape_plot_woCNAs              <- args$output_path_landscape_plot_woCNAs
output_path_CNA_comparison_plot	               <- args$output_path_CNA_comparison_plot
output_path_CNA_comparison_plot_woCNAs	       <- args$output_path_CNA_comparison_plot_woCNAs

`%>%` <- magrittr::`%>%`

metadata <- read.table(input_path_metadata, header = T, sep = "\t")

counts        <- read.table(input_path_counts, header = T, sep = "\t")
counts_woCNAs <- read.table(input_path_counts_woCNAs, header = T, sep = "\t")
#counts_sva <- read.table(input_path_counts_sva, header = T, sep = "\t")

library_sizes        <- read.table(input_path_library_sizes, header = T, sep = "\t")
library_sizes_woCNAs <- read.table(input_path_library_sizes_woCNAs, header = T, sep = "\t")
#library_sizes <- library_sizes[match(colnames(counts), library_sizes$SRA_id),]

load(input_path_model_matrix); mm_wCNAs <- mm
load(input_path_model_matrix_woCNAs); mm_woCNAs <- mm

ampl <- read.table(input_path_ampl, header = T, sep = "\t")
del <- read.table(input_path_del, header = T, sep = "\t")

out.list <- do_DA_analysis(counts = counts, library_sizes = library_sizes, mm = mm_wCNAs)
out.list.woCNAs <- do_DA_analysis(counts = counts_woCNAs, library_sizes = library_sizes_woCNAs, mm = mm_woCNAs)

write.table(x = out.list$topTags, file = output_path_DARs, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(x = out.list.woCNAs$topTags, file = output_path_DARs_woCNAs, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

ggplot2::ggsave(output_path_diagnostics_plot, out.list$p_diagnostics, width = 8, height = 6)
ggplot2::ggsave(output_path_volcano_plot, out.list$p_volcano, width = 6, height = 5)
ggplot2::ggsave(output_path_landscape_plot, out.list$p_landscape, width = 10, height = 4)
ggplot2::ggsave(output_path_CNA_comparison_plot, out.list$p_CNA_comparison, width = 10, height = 6)

ggplot2::ggsave(output_path_diagnostics_plot_woCNAs, out.list.woCNAs$p_diagnostics, width = 8, height = 6)
ggplot2::ggsave(output_path_volcano_plot_woCNAs, out.list.woCNAs$p_volcano, width = 6, height = 5)
ggplot2::ggsave(output_path_landscape_plot_woCNAs, out.list.woCNAs$p_landscape, width = 10, height = 4)
ggplot2::ggsave(output_path_CNA_comparison_plot_woCNAs, out.list.woCNAs$p_CNA_comparison, width = 10, height = 6)
