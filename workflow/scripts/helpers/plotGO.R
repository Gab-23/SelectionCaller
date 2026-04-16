plotGO <- function(df){

p <- ggplot2::ggplot(df, ggplot2::aes(x = Enrichment, y = Term)) +
  ggplot2::geom_point(ggplot2::aes(size = sizeVar, color = signif)) +
  ggplot2::scale_color_manual(values = c("Significant" = "red", "Not significant" = "grey60")) +
  ggplot2::scale_size_continuous(name = "-log10(FDR)") +
  ggplot2::theme_bw() +
  ggplot2::labs(
    x = "Enrichment",
    y = "GO Term",
    color = "FDR < 0.05"
  ) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(face = "bold"),
    axis.text.y = ggplot2::element_text(size = 8)
  )

return(p)

}