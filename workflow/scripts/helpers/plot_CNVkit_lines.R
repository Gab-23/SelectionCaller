plot_CNVkit_lines <- function(df){

df.out <- df %>%
  dplyr::mutate(
    chromosome = sub("_.*", "", binID),
    bin = as.numeric(sub(".*_", "", binID)),
    chromosome = factor(chromosome, levels = chr_levels)
  )    

df.out$binID <- NULL

df.out.long <- df.out %>%
  dplyr::pivot_longer(
    cols = -c(bin, chromosome),   # keep these, reshape the rest
    names_to = "sample",
    values_to = "scaling_factors"
  )

p <- ggplot2::ggplot(df.out.long,
                ggplot2::aes(x = bin, y = scaling_factors,
                             group = sample, color = sample)) +
  ggplot2::geom_line(alpha = 0.4) +
  ggplot2::stat_summary(
    ggplot2::aes(group = 1),
    fun = mean,
    geom = "line",
    color = "black",
    linewidth = 1.2
  ) +                               
  ggplot2::facet_wrap(~ chromosome, scales = "free_x", nrow = 1) +
  ggplot2::labs(                                                  
    x = "Genomic bin",                                            
    y = "Scaling factor",
    title = ""
  ) +                                         
  ggplot2::theme_bw() +                                         
  ggplot2::theme(      
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold"),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
  ) + 
  ggplot2::ylim(0, 3)

return(p)

}
