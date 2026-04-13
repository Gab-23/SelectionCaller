plot_landscape <- function(df, what_to_plot, y_axis = "Generic ES", title = "Generic title"){
  
  p <- ggplot2::ggplot(df,
  ggplot2::aes(x = bin, y = .data[[what_to_plot]], group = chromosome)) +
  ggplot2::geom_line(alpha = 0.8) +
  ggplot2::facet_wrap(~ chromosome, scales = "free_x", nrow = 1) +
  ggplot2::labs(
    x = "Genomic bin",
    y = y_axis,
    title = title
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold"),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
  )

}
