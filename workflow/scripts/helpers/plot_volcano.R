plot_volcano <- function(topTags, top_n = 10){
	
	topDARs <- topTags[order(topTags$FDR, -abs(topTags$logFC)), ][1:top_n, ]

	p_volcano <- ggplot2::ggplot(topTags, ggplot2::aes(x = logFC, y = -log10(FDR))) +
	  ggplot2::geom_point(ggplot2::aes(color = significant), size = 0.7) +
	  ggplot2::scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "red")) +
	  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
	  ggplot2::geom_text(
		data = topDARs,
		ggplot2::aes(x = logFC, y = -log10(FDR), label = binID),
		color = "blue",
		size = 3,
		vjust = -0.5
	  ) +
	  ggplot2::ggtitle("Volcano plot (top 10 DARs labeled)") +
	  ggplot2::xlab("logFC") +
	  ggplot2::ylab("-log10(FDR)") +
	  ggplot2::theme_bw() +
	  ggplot2::theme(legend.title = ggplot2::element_blank())

}
