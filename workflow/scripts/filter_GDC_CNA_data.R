options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path                       <- args$input_path
output_path_filtered             <- args$output_path_filtered
output_path_plot_probes          <- args$output_path_plot_probes
output_path_plot_segment         <- args$output_path_plot_segment
min_Num_Probes                   <- as.numeric(args$min_Num_Probes)
min_abs_Segment_Mean             <- as.numeric(args$min_abs_Segment_Mean)

segmentation_file <- read.delim(input_path, header = F)
# NOTE (29/10/2026): l'aggregato ha ora una sesta colonna, l'identificativo del campione (aggregate_GDC_CNA_data), che serve
# a contare i pazienti distinti per bin; viene mantenuta nel file filtrato. Accetta anche il vecchio formato a cinque colonne.
if (ncol(segmentation_file) == 6) {
    colnames(segmentation_file) <- c("chr", "start", "end", "num_probes", "segment_mean", "sample")
} else {
    colnames(segmentation_file) <- c("chr", "start", "end", "num_probes", "segment_mean")
}

segmentation_file.f <- segmentation_file[(segmentation_file$num_probes > min_Num_Probes) & 
                                         (abs(segmentation_file$segment_mean) > min_abs_Segment_Mean),]


boxplot_data <- rbind(
  data.frame(
    num_probes = segmentation_file$num_probes,
    filtering = "Pre filtering"
  ),
  data.frame(
    num_probes = segmentation_file.f$num_probes,
    filtering = "Post filtering"
  )
)

boxplot_data$filtering <- factor(
  boxplot_data$filtering,
  levels = c("Pre filtering", "Post filtering")
)

hist_data <- rbind(
  data.frame(
    segment_mean = segmentation_file$segment_mean,
    filtering = "Pre filtering"
  ),
  data.frame(
    segment_mean = segmentation_file.f$segment_mean,
    filtering = "Post filtering"
  )
)

hist_data$filtering <- factor(
  hist_data$filtering,
  levels = c("Pre filtering", "Post filtering")
)

box <- ggplot2::ggplot(
  boxplot_data,
  ggplot2::aes(
    x = filtering,
    y = num_probes
  )
) +
  ggplot2::geom_boxplot() +
  ggplot2::labs(
    x = NULL,
    y = "num_probes",
    title = "Number of probes per segment"
  ) +
  ggplot2::theme_minimal()

ggplot2::ggsave(
  filename = output_path_plot_probes,
  plot = box,
  device = "pdf"
)

hist <- ggplot2::ggplot(
  hist_data,
  ggplot2::aes(x = segment_mean)
) +
  ggplot2::geom_histogram(bins = 50) +
  ggplot2::facet_wrap(~ filtering, ncol = 2) +
  ggplot2::labs(
    x = "segment_mean",
    y = "Count",
    title = "Distribution of segment means"
  ) +
  ggplot2::theme_minimal()

ggplot2::ggsave(
  filename = output_path_plot_segment,
  plot = hist,
  device = "pdf"
)

write.table(x = segmentation_file.f,                                              # write output
                file = output_path_filtered, 
                quote = F, 
                sep = "\t", 
                row.names = F, 
                col.names = F)