source("workflow/scripts/helpers/edg2ATAC.R")
source("workflow/scripts/helpers/plot_diagnostics.R")
source("workflow/scripts/helpers/plot_volcano.R")
source("workflow/scripts/helpers/plot_landscape_FC.R")
source("workflow/scripts/helpers/compute_FC_manually.R")
source("workflow/scripts/helpers/process_topTags.R")
source("workflow/scripts/helpers/plot_comparison_FC.R")
source("workflow/scripts/helpers/compare_acc_CNAs.R")

library(edgeR)
library(sva)

options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_counts                     <- args$input_path_counts
input_path_counts_sva	              <- args$input_path_counts_sva
input_path_metadata                   <- args$input_path_metadata
input_path_library_sizes              <- args$input_path_library_sizes
input_path_model_matrix               <- args$input_path_model_matrix
input_path_ampl			      <- args$input_path_ampl
input_path_del                        <- args$input_path_del

output_path_DARs		      <- args$output_path_DARs
output_path_manual_DARs		      <- args$output_path_manual_DARs
output_path_diagnostics_plot	      <- args$output_path_diagnostics_plot
output_path_volcano_plot              <- args$output_path_volcano_plot
output_path_landscape_plot            <- args$output_path_landscape_plot
output_path_manual_landscape_plot     <- args$output_path_manual_landscape_plot
output_path_comparison_plot	      <- args$output_path_comparison_plot
output_path_CNA_comparison_plot	      <- args$output_path_CNA_comparison_plot

`%>%` <- magrittr::`%>%`

counts <- read.table(input_path_counts, header = T, sep = "\t")
counts_sva <- read.table(input_path_counts_sva, header = T, sep = "\t")

metadata <- read.table(input_path_metadata, header = T, sep = "\t")

library_sizes <- read.table(input_path_library_sizes, header = T, sep = "\t")
#library_sizes <- library_sizes[match(colnames(counts), library_sizes$SRA_id),]

load(input_path_model_matrix)

ampl <- read.table(input_path_ampl, header = T, sep = "\t")
del <- read.table(input_path_del, header = T, sep = "\t")

rownames(counts) <- counts$binID
counts$binID <- NULL

rownames(counts_sva) <- counts_sva$binID
counts_sva$binID <- NULL

libsize <- library_sizes$library_size
names(libsize) <- library_sizes$SRA_id

names(libsize) == colnames(counts)

topTags <- edg2ATAC(counts, mm, libsize, "conditiontumor")
topTags$binID <- rownames(topTags)
rownames(topTags) <- NULL

manual_topTags <- compute_FC_manually(counts_sva, metadata)
colnames(manual_topTags) <- c("binID", "NAT", "tumor", "FC", "logFC")
manual_topTags$PValue <- rep(1, nrow(manual_topTags))
manual_topTags$FDR <- rep(1, nrow(manual_topTags))

topTags <- process_topTags(topTags)
manual_topTags <- process_topTags(manual_topTags)
  
p_diagnostics <- plot_diagnostics(topTags)
p_volcano <- plot_volcano(topTags)
p_landscape <- plot_landscape_FC(topTags)
p_landscape_manual <- plot_landscape_FC(manual_topTags)

p_comparison <- plot_comparison_FC(topTags, manual_topTags)
p_CNA_comparison <- compare_acc_CNAs(topTags, ampl, del)

ggplot2::ggsave(output_path_diagnostics_plot, p_diagnostics, width = 8, height = 6)
ggplot2::ggsave(output_path_volcano_plot, p_volcano, width = 6, height = 5)
ggplot2::ggsave(output_path_landscape_plot, p_landscape, width = 10, height = 4)
ggplot2::ggsave(output_path_manual_landscape_plot, p_landscape_manual, width = 10, height = 4)
ggplot2::ggsave(output_path_comparison_plot, p_comparison, width = 6, height = 6)
ggplot2::ggsave(output_path_CNA_comparison_plot, p_CNA_comparison, width = 10, height = 6)

write.table(x = topTags, 
	    file = output_path_DARs, 
	    sep = "\t",
            col.names = TRUE, 
            row.names = FALSE, 
            quote = FALSE)

write.table(x = manual_topTags,
	    file = output_path_manual_DARs,
	    sep = "\t",
	    col.names = TRUE,
	    row.names = FALSE,
	    quote = FALSE)

