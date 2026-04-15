plot_landscape_FC <- function(topTags){
	
	chr_levels <- sort(unique(topTags$chr), na.last = TRUE)
	chr_levels <- chr_levels[order(as.numeric(gsub("chr", "", chr_levels)))]
	topTags$chr <- factor(topTags$chr, levels = chr_levels)

	chr_lengths <- tapply(topTags$start, topTags$chr, max)
	chr_offsets <- cumsum(c(0, head(chr_lengths, -1)))
	names(chr_offsets) <- names(chr_lengths)
	
	topTags$cum_pos <- topTags$start + chr_offsets[topTags$chr]

	chr_boundaries <- data.frame(
	  chr = names(chr_lengths),
	  start = chr_offsets,
	  end = chr_offsets + chr_lengths
	)

	chr_centers <- data.frame(
	  chr = names(chr_lengths),
	  center = chr_offsets + chr_lengths / 2
	)

	p_fc <- ggplot2::ggplot(topTags, ggplot2::aes(x = cum_pos, y = logFC)) +
	  ggplot2::geom_point(ggplot2::aes(color = significant), size = 0.5) +
	  ggplot2::scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "red")) +
	  ggplot2::geom_vline(data = chr_boundaries, ggplot2::aes(xintercept = start),
		                  linetype = "dashed", size = 0.3) +
	  ggplot2::scale_x_continuous(breaks = chr_centers$center, labels = chr_centers$chr) +
	  ggplot2::ggtitle("Genomic logFC (FDR < 0.05 highlighted)") +
	  ggplot2::xlab("Chromosome") +
	  ggplot2::ylab("logFC") +
	  ggplot2::theme_bw() +
	  ggplot2::theme(
		legend.title = ggplot2::element_blank(),
		axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5)
	  )
}

