plot_barplot <- function(df, col1, col2, y_axis = "Generic ES"){

p.barplot_samples <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[col1]], y = .data[[col2]], fill = .data[[col1]])) +
  ggplot2::geom_bar(stat = "identity", color = "black") +
  ggplot2::labs(
    title = "",
    x = "",
    y = y_axis
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "none")

}