source("workflow/scripts/helpers/plot_landscape.R")

options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_GC                       <- args$input_path_ampl
input_path_models_DepMap            <- args$input_path_models_DepMap
input_path_manifest_GDC             <- args$input_path_manifest_GDC
input_path_samples_COSMIC           <- args$input_path_samples_COSMIC
output_path_GC_landscape            <- args$output_path_GC_landscape
output_path_GC_distribution         <- args$output_path_GC_distribution
output_path_samples_barplot         <- args$output_path_samples_barplot

`%>%` <- magrittr::`%>%`

GC <- read.table(input_path_GC, header = T, sep = "\t")
colnames(GC) <- c("chr", "start", "end", "binID", "GC.content")

DepMap <- readLines(input_path_models_DepMap)
COSMIC <- readLines(input_path_samples_COSMIC)
GDC <- read.table(input_path_manifest_GDC, header = T, sep = "\t")

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

p.GC <- plot_landscape(GC.out, "GC_content", "GC content", "GC distribution")

n_samples_DepMap <- length(DepMap)
n_samples_COSMIC <- length(COSMIC)
n_samples_GDC <- nrow(GDC)

samples_df <- data.frame(
  Dataset = c("DepMap", "COSMIC", "GDC"),
  Samples = c(n_samples_DepMap, n_samples_COSMIC, n_samples_GDC)
)

p.barplot <- ggplot2::ggplot(samples_df, ggplot2::aes(x = Dataset, y = Samples, fill = Dataset)) +
  ggplot2::geom_bar(stat = "identity", color = "black") +
  ggplot2::labs(
    title = "",
    x = "",
    y = "# Samples"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "none")

ggplot2::ggsave(output_path_GC_landscape, p.GC, width = 14, height = 4)
ggplot2::ggsave(output_path_GC_distribution, p.hist, width = 6, height = 4)
ggplot2::ggsave(output_path_samples_barplot, p.barplot, width = 6, height = 4)





