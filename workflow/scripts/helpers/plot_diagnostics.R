plot_diagnostics <- function(topTags){
  
  plot_df <- data.frame(
  value = c(topTags$logCPM, topTags$logFC, topTags$PValue, topTags$FDR),
  variable = factor(
    rep(c("logCPM", "logFC", "PValue", "FDR"), each = nrow(topTags)),
    levels = c("logCPM", "logFC", "PValue", "FDR")
  )
)

p_all <- ggplot2::ggplot(plot_df, ggplot2::aes(x = value)) +
  ggplot2::geom_histogram(bins = 50) +
  ggplot2::facet_wrap(~variable, scales = "free", ncol = 2) +
  ggplot2::ggtitle("Distributions: logCPM, logFC, PValue, FDR") +
  ggplot2::theme_bw()
  
}

