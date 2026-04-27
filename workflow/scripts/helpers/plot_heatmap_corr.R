plot_heatmap_corr <- function(cor.df){

    p <- ggplot2::ggplot(cor.df, ggplot2::aes(x = Var1, y = Var2, fill = Correlation)) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    axis.title = ggplot2::element_blank()
  )
  
}