source("workflow/scripts/helpers/plot_landscape.R")
source("workflow/scripts/helpers/plot_barplot.R")

options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_GC                                 <- args$input_path_GC
input_path_models_DepMap                      <- args$input_path_models_DepMap
input_path_manifest_GDC                       <- args$input_path_manifest_GDC
input_path_samples_COSMIC                     <- args$input_path_samples_COSMIC
input_path_OGs                                <- args$input_path_OGs
input_path_TSGs                               <- args$input_path_TSGs
input_path_essential                          <- args$input_path_essential
input_path_coding_mutations                   <- args$input_path_coding_mutations
input_path_noncoding_mutations                <- args$input_path_noncoding_mutations
output_path_GC_landscape                      <- args$output_path_GC_landscape
output_path_GC_distribution                   <- args$output_path_GC_distribution
output_path_samples_barplot                   <- args$output_path_samples_barplot
output_path_genes_barplot                     <- args$output_path_genes_barplot
output_path_num_mutations_barplot             <- args$output_path_num_mutations_barplot
output_path_coding_mut_descriptions_barplot   <- args$output_path_coding_mut_descriptions_barplot

`%>%` <- magrittr::`%>%`

GC <- read.table(input_path_GC, header = F, sep = "\t")
colnames(GC) <- c("chr", "start", "end", "binID", "GC.content")

DepMap <- readLines(input_path_models_DepMap)
COSMIC <- readLines(input_path_samples_COSMIC)
GDC <- read.table(input_path_manifest_GDC, header = T, sep = "\t")

coding <- read.table(input_path_coding_mutations, header = F, sep = "\t")
colnames(coding) <- c("chr", "start", "end", "mut_id", "description")

noncoding <- read.table(input_path_noncoding_mutations, header = F, sep = "\t")
colnames(noncoding) <- c("chr", "start", "end", "mut_id")

OGs <- read.table(input_path_OGs, header = F, sep = "\t")
colnames(OGs) <- c("chr", "start", "end", "TxID", "gene_name", "gene_class")

TSGs <- read.table(input_path_TSGs, header = F, sep = "\t")
colnames(TSGs) <- c("chr", "start", "end", "TxID", "gene_name", "gene_class")

essential <- read.table(input_path_essential, header = F, sep = "\t")
colnames(essential) <- c("chr", "start", "end", "gene_name")

p.hist <- ggplot2::ggplot(GC, ggplot2::aes(x = GC.content)) +
  ggplot2::geom_histogram(
    bins = 50, 
    fill = "steelblue", 
    color = "black"
  ) +
  ggplot2::labs(
    title = "GC Content Histogram",
    x = "GC content",
    y = "Frequency"
  ) +
  ggplot2::theme_minimal()

chr_levels <- 1:22

GC.out <- GC %>%
  dplyr::select(binID, GC.content) %>%
  dplyr::mutate(
    chromosome = sub("_.*", "", binID),
    bin = as.numeric(sub(".*_", "", binID)),
    chromosome = factor(chromosome, levels = chr_levels)
  )

p.GC <- plot_landscape(GC.out, "GC.content", "GC content", "GC distribution")

n_samples_DepMap <- length(DepMap)
n_samples_COSMIC <- length(COSMIC)
n_samples_GDC <- nrow(GDC)

samples_df <- data.frame(
  Dataset = c("DepMap", "COSMIC", "GDC"),
  Samples = c(n_samples_DepMap, n_samples_COSMIC, n_samples_GDC)
)

p.barplot_samples <- plot_barplot(samples_df, "Dataset", "Samples", "# Samples")

n_OGs <- nrow(OGs)
n_TSGs <- nrow(TSGs)
n_essential <- nrow(essential)

genes_df <- data.frame(
  Category = c("OGs", "TSGs", "essential"),
  Genes = c(n_OGs, n_TSGs, n_essential)
)

p.barplot_genes <- plot_barplot(genes_df, "Category", "Genes", "# Genes")

n_coding_mutations <- nrow(coding)
n_noncoding_mutations <- nrow(noncoding)

n_mut_df <- data.frame(
  Category = c("coding", "noncoding"),
  Mutations = c(n_coding_mutations, n_noncoding_mutations)
)

p.barplot_n_mut <- plot_barplot(n_mut_df, "Category", "Mutations", "# Mutations")

n_mut_by_description <- coding %>%
  tidyr::separate_rows(description, sep = ",") %>%
  dplyr::group_by(description) %>%
  dplyr::summarise(Mutations = dplyr::n(), .groups = "drop") %>%
  dplyr::rename(Category = description)

p.barplot_n_mut_by_description <- plot_barplot(n_mut_by_description, "Category", "Mutations", "# Mutations") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

ggplot2::ggsave(output_path_GC_landscape, p.GC, width = 14, height = 4)
ggplot2::ggsave(output_path_GC_distribution, p.hist, width = 6, height = 4)
ggplot2::ggsave(output_path_samples_barplot, p.barplot_samples, width = 6, height = 4)
ggplot2::ggsave(output_path_genes_barplot, p.barplot_genes, width = 6, height = 4)
ggplot2::ggsave(output_path_num_mutations_barplot, p.barplot_n_mut, width = 6, height = 4)
ggplot2::ggsave(output_path_coding_mut_descriptions_barplot, p.barplot_n_mut_by_description, width = 6, height = 4)