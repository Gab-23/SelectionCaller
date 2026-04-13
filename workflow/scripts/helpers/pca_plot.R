pca_plot <- function(pca, pca_df, hue, label = F, shape = NULL) {

  if (is.null(shape)) {
    p <- ggplot2::ggplot(
      pca_df,
      ggplot2::aes(x = PC1, y = PC2, color = .data[[hue]])
      ) +
      ggplot2::geom_point(size = 3)
  } else {
    p <- ggplot2::ggplot(
      pca_df,
      ggplot2::aes(
        x = PC1,
        y = PC2,
        color = .data[[hue]],
        shape = .data[[shape]]
       )
      ) +
      ggplot2::geom_point(size = 3)
  }

  if (label){
    p <- p + ggplot2::geom_text(ggplot2::aes(label = SRA_id), vjust = -0.5, inherit.aes = TRUE, show.legend = FALSE)
  }

  p + 
  ggplot2::theme_bw() +
  ggplot2::labs(
    title = paste0(""),
    x = paste0("PC1 (", round(summary(pca)$importance[2,1] * 100, 2), "%)"),
    y = paste0("PC2 (", round(summary(pca)$importance[2,2] * 100, 2), "%)")
  )
}



