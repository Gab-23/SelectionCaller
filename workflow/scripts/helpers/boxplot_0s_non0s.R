boxplots_0s_non0s <- function(df, pattern){

  cols_to_plot <- grep(pattern = pattern, x = colnames(df), value = T) # "length_(non)?coding"

  df <- df[,c("logFC", cols_to_plot)]
  
  df_long <- df %>%
    tidyr::pivot_longer(
      cols = -logFC,
      names_to = "variable",
      values_to = "value"
    ) %>%
    dplyr::mutate(group = ifelse(value == 0, "0", "non-0"))

  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = group, y = logFC, fill = group)) +
    ggplot2::geom_boxplot() +
    ggplot2::facet_wrap(~variable, scales = "free") +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = "Value group",
      y = "logFC",
      title = "logFC distribution: 0 vs non-0 across variables"
    )
       
  return(p)
}