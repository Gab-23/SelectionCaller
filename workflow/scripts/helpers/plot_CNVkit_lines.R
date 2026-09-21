plot_CNVkit_lines <- function(df, metadata){

chr_levels <- 1:22

df.out <- df %>%
  dplyr::mutate(
    chromosome = sub("_.*", "", binID),
    bin = as.numeric(sub(".*_", "", binID)),
    chromosome = factor(chromosome, levels = chr_levels)
  )    

df.out$binID <- NULL

df.out.long <- df.out %>%
  tidyr::pivot_longer(
    cols = -c(bin, chromosome),
    names_to = "sample",
    values_to = "scaling_factors"
  )

df.out.long <- dplyr::left_join(
  df.out.long,
  metadata[, c("SRA_id", "condition")],
  by = c("sample" = "SRA_id")
)

if (any(is.na(df.out.long$condition))) {
  missing_samples <- unique(df.out.long$sample[is.na(df.out.long$condition)])
  warning(paste0("plot_CNVkit_lines: no condition found in metadata for sample(s), excluded from plot: ",
                 paste(missing_samples, collapse = ", ")))
  df.out.long <- df.out.long[!is.na(df.out.long$condition), ]
}

p <- ggplot2::ggplot(df.out.long,
                ggplot2::aes(x = bin, y = scaling_factors,
                             group = sample, color = condition)) +
  ggplot2::stat_summary(
    ggplot2::aes(group = condition, color = condition),
    fun = mean,
    geom = "line",
    linewidth = 0.7
  ) +
  ggplot2::facet_wrap(~ chromosome, scales = "free_x", nrow = 1) +
  ggplot2::labs(
    x = "Genomic bin",
    y = "Scaling factor",
    color = "Condition",
    title = ""
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold"),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
  ) +
  ggplot2::ylim(0, 3)

return(p)

}