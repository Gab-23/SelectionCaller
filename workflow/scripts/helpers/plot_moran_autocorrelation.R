plot_moran_autocorrelation <- function(df, K){

  plot_autocorrelation <- ggplot2::ggplot(                         
			        df,
			        ggplot2::aes(x = x, y = y)
				) +
			        ggplot2::geom_point(size = 1.5, alpha = 0.7) +
			        ggplot2::geom_smooth(method = "lm", se = FALSE, color = "black") +
				ggplot2::geom_vline(xintercept = mean(df$x, na.rm = TRUE), linetype = "dashed") +
				ggplot2::geom_hline(yintercept = mean(df$y, na.rm = TRUE), linetype = "dashed") +
				ggplot2::geom_text(ggplot2::aes(label = annot), vjust = -0.5, size = 3) +
				ggplot2::labs(
				      x = paste0("I round turnout [", K, "]"),
			              y = "Lagged turnout",
				      title = "Moran Scatterplot"
			              ) +
				ggplot2::theme_bw()

  return(plot_autocorrelation)

}
