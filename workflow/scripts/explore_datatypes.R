source("workflow/scripts/helpers/plot_landscape.R")

options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_GC                       <- args$input_path_ampl
output_path_GC_landscape            <- args$output_path_GC_landscape
output_path_GC_distribution         <- args$output_path_GC_distribution

`%>%` <- magrittr::`%>%`

GC <- read.table(input_path_GC, header = T, sep = "\t")
colnames(GC) <- c("chr", "start", "end", "binID", "GC.content")

p.hist <- ggplot2::ggplot(GC, ggplot2::aes(x = GC_content)) +
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

ggplot2::ggsave(output_path_GC_landscape, p.GC, width = 14, height = 4)
ggplot2::ggsave(output_path_GC_distribution, p.hist, width = 6, height = 4)
