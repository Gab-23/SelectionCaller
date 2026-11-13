source("workflow/scripts/helpers/edg2ATAC.R")
source("workflow/scripts/helpers/plot_diagnostics.R")
source("workflow/scripts/helpers/plot_volcano.R")
source("workflow/scripts/helpers/plot_landscape_FC.R")
#source("workflow/scripts/helpers/compute_FC_manually.R")
source("workflow/scripts/helpers/process_topTags.R")
source("workflow/scripts/helpers/plot_comparison_FC.R")
source("workflow/scripts/helpers/compare_acc_CNAs.R")
source("workflow/scripts/helpers/do_DA_analysis.R")
source("workflow/scripts/helpers/plot_GC_bias.R")                                 # NEW: GC-bias diagnostic
source("workflow/scripts/helpers/plot_GC_curves_per_sample.R")                    # NEW: per-sample GC curves
source("workflow/scripts/helpers/compute_GC_offset.R")                            # NEW: optional per-sample GC offset

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
input_path_GC                                  <- args$input_path_GC              # NEW: binID, GC_content (process_GC_content)

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
output_path_GC_bias_plot                       <- args$output_path_GC_bias_plot          # NEW
output_path_GC_bias_plot_woCNAs                <- args$output_path_GC_bias_plot_woCNAs   # NEW
output_path_GC_curves_plot                     <- args$output_path_GC_curves_plot        # NEW
output_path_GC_curves_plot_woCNAs              <- args$output_path_GC_curves_plot_woCNAs # NEW

# NEW: GC correction parameters (config: GC_correction_atacseqs, GC_span_atacseqs, GC_n_fit_atacseqs)
GC_correction                                  <- args$GC_correction                      # "none" | "loess"
GC_span                                        <- as.numeric(args$GC_span)
GC_n_fit                                       <- as.integer(args$GC_n_fit)

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

gc <- read.table(input_path_GC, header = T, sep = "\t")                          # NEW

out.list <- do_DA_analysis(counts = counts, library_sizes = library_sizes, mm = mm_wCNAs,
                           gc = gc, gc_correction = GC_correction, gc_span = GC_span, gc_n_fit = GC_n_fit)          # NEW: GC args
out.list.woCNAs <- do_DA_analysis(counts = counts_woCNAs, library_sizes = library_sizes_woCNAs, mm = mm_woCNAs,
                                  gc = gc, gc_correction = GC_correction, gc_span = GC_span, gc_n_fit = GC_n_fit)   # NEW: GC args

GC_tag <- if (GC_correction == "none") "no GC correction" else sprintf("GC correction: %s, span %g", GC_correction, GC_span)   # NEW

# NEW: GC-bias diagnostic on the tested bins (logFC and logCPM by GC-quantile bin)
# (topTags are AFTER correction, so this plot is also the check that the correction worked)
p_GC_bias        <- plot_GC_bias(out.list$topTags, gc, n_bins = 20, title_suffix = paste0(" (wCNAs, ", GC_tag, ")"))
p_GC_bias_woCNAs <- plot_GC_bias(out.list.woCNAs$topTags, gc, n_bins = 20, title_suffix = paste0(" (woCNAs, ", GC_tag, ")"))

# NEW: per-sample GC curves on the same counts/library sizes that enter the DA analysis
# (colour = condition, columns = study_id; subtitle: slope variance explained by study vs condition)
# page 1: before correction; page 2 (only if a correction is applied): after subtracting the GC offset
GC_curves_pages <- function(cnt, libs, offset, label){
  pages <- list(plot_GC_curves_per_sample(cnt, libs, metadata, gc, n_bins = 20, title_suffix = paste0(" (", label, ", before GC correction)")))
  if (!is.null(offset)) pages[[2]] <- plot_GC_curves_per_sample(cnt, libs, metadata, gc, n_bins = 20,
                                                                 title_suffix = paste0(" (", label, ", after ", GC_tag, ")"), gc_offset = offset)
  pages
}
p_GC_curves        <- GC_curves_pages(counts, library_sizes, out.list$gc_offset, "wCNAs")
p_GC_curves_woCNAs <- GC_curves_pages(counts_woCNAs, library_sizes_woCNAs, out.list.woCNAs$gc_offset, "woCNAs")
save_pages <- function(pages, path){ grDevices::pdf(path, width = 12, height = 8); for (p in pages) print(p); invisible(grDevices::dev.off()) }

write.table(x = out.list$topTags, file = output_path_DARs, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(x = out.list.woCNAs$topTags, file = output_path_DARs_woCNAs, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

ggplot2::ggsave(output_path_diagnostics_plot, out.list$p_diagnostics, width = 8, height = 6)
ggplot2::ggsave(output_path_volcano_plot, out.list$p_volcano, width = 6, height = 5)
ggplot2::ggsave(output_path_landscape_plot, out.list$p_landscape, width = 10, height = 4)
ggplot2::ggsave(output_path_CNA_comparison_plot, out.list$p_CNA_comparison, width = 10, height = 6)
ggplot2::ggsave(output_path_GC_bias_plot, p_GC_bias, width = 10, height = 8)                    # NEW
save_pages(p_GC_curves, output_path_GC_curves_plot)                                             # NEW

ggplot2::ggsave(output_path_diagnostics_plot_woCNAs, out.list.woCNAs$p_diagnostics, width = 8, height = 6)
ggplot2::ggsave(output_path_volcano_plot_woCNAs, out.list.woCNAs$p_volcano, width = 6, height = 5)
ggplot2::ggsave(output_path_landscape_plot_woCNAs, out.list.woCNAs$p_landscape, width = 10, height = 4)
ggplot2::ggsave(output_path_CNA_comparison_plot_woCNAs, out.list.woCNAs$p_CNA_comparison, width = 10, height = 6)
ggplot2::ggsave(output_path_GC_bias_plot_woCNAs, p_GC_bias_woCNAs, width = 10, height = 8)      # NEW
save_pages(p_GC_curves_woCNAs, output_path_GC_curves_plot_woCNAs)                               # NEW