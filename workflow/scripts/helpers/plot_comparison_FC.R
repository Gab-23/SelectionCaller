plot_comparison_FC <- function(topTags, manual_topTags){

  df <- merge(x = topTags[,c("binID", "logFC")], 
  			  y = manual_topTags[,c("binID", "logFC")],
  			  by = "binID",
  			  suffixes = c("","_manual"))

  pearson_cor  <- stats::cor(df$logFC, df$logFC_manual, method = "pearson")
  spearman_cor <- stats::cor(df$logFC, df$logFC_manual, method = "spearman")

  label_text <- paste0(
    "Pearson = ", round(pearson_cor, 3),
    "\nSpearman = ", round(spearman_cor, 3)
  )

  x_pos <- min(df$logFC, na.rm = TRUE)
  y_pos <- max(df$logFC_manual, na.rm = TRUE)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = logFC, y = logFC_manual)) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
    ggplot2::annotate(
      "text",
      x = x_pos,
      y = y_pos,
      label = label_text,
      hjust = 0,
      vjust = 1,
      size = 4
    ) +
    ggplot2::labs(
      x = "logFC (topTags)",
      y = "logFC (manual_topTags)",
      title = "logFC Comparison"
    ) +
    ggplot2::theme_bw()

  return(p)
}
