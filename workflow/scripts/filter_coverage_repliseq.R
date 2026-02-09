options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path                       <- args$input_path
output_path_filtered             <- args$output_path_filtered
output_path_plot                 <- args$output_path_plot
cov_thr                          <- as.numeric(args$cov_thr)

coverage <- read.table(file = input_path, header = F, sep = "\t")
colnames(coverage) <- c("chr", "start", "end", "binID", "coverage")

coverage.f <- coverage[coverage$coverage > cov_thr,] # based on G1 50

qs <- quantile(
  coverage.f$coverage,
  probs = c(0.01, 0.025, 0.25, 0.5, 0.75, 0.975, 0.99),
  na.rm = TRUE
)

hist <- ggplot2::ggplot(coverage.f, ggplot2::aes(x = coverage)) +
  ggplot2::geom_histogram(
    binwidth = 1,
    fill = "skyblue",
    color = "black",
    alpha = 0.7
  ) +
  ggplot2::geom_vline(
    xintercept = qs,
    linetype = "dashed",
    color = "red",
    linewidth = 0.8
  ) +
  ggplot2::labs(
    title = "Raw Coverage",
    x = "Read Count",
    y = "Frequency"
  ) +
  ggplot2::scale_x_continuous(n.breaks = 30) +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
  )  

ggplot2::ggsave(
  filename = output_path_plot,
  plot = hist,
  width = 8,
  height = 6,
  dpi = 300
)

coverage.ff <- coverage.f[(coverage.f$coverage > qs[2]) &                     # filter on bilateral 2.5 %-ile Repli-scan style
                          (coverage.f$coverage < qs[6]),]

write.table(x = coverage.ff,                                                  # write output
            file = output_path_filtered, 
            quote = F, 
            sep = "\t", 
            row.names = F, 
            col.names = F)