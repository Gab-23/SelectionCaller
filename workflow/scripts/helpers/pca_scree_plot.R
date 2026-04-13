pca_scree_plot <- function(variance_vector){

df <- data.frame(
  PC = factor(paste0("PC", seq_along(variance_vector)), levels = paste0("PC", seq_along(variance_vector))),
  Variance = variance_vector
)

p <- ggplot2::ggplot(df, ggplot2::aes(x = as.numeric(PC), y = Variance)) +
  ggplot2::geom_line() +
  ggplot2::geom_point() +
  ggplot2::theme_bw() +
  ggplot2::labs(
    title = "Scree Plot",
    x = "Principal Component",
    y = "Proportion of Variance Explained"
  )

}
