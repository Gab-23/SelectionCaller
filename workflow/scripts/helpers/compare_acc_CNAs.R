compare_acc_CNAs <- function(topTags, ampl, del){

  FC_w_ampl <- merge(topTags[,c("binID", "chr", "logFC")], ampl, by = "binID")
  FC_w_del  <- merge(topTags[,c("binID", "chr", "logFC")], del, by = "binID")

  FC_w_ampl$type <- "ampl"
  FC_w_del$type  <- "del"

  colnames(FC_w_ampl) <- c("binID", "chr", "logFC", "frequency", "type")
  colnames(FC_w_del)  <- c("binID", "chr", "logFC", "frequency", "type")

  df <- base::rbind(FC_w_ampl, FC_w_del)

  cor_df <- do.call(rbind, lapply(split(df, df$type), function(subdf){

    pearson  <- stats::cor(subdf$frequency, subdf$logFC, method = "pearson", use = "complete.obs")
    spearman <- stats::cor(subdf$frequency, subdf$logFC, method = "spearman", use = "complete.obs")

    data.frame(
      type = unique(subdf$type),
      label = paste0(
        "Pearson = ", round(pearson, 3),
        "\nSpearman = ", round(spearman, 3)
      ),
      x = min(subdf$frequency, na.rm = TRUE),
      y = max(subdf$logFC, na.rm = TRUE)
    )
  }))

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = frequency, y = logFC, color = chr)
  ) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "black") +
    ggplot2::facet_wrap(~type, scales = "free_x") +
    ggplot2::geom_text(
      data = cor_df,
      ggplot2::aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 4
    ) +
    ggplot2::labs(
      x = "Frequency (ampl / del)",
      y = "logFC",
      title = "logFC vs Frequency (Amplifications vs Deletions)",
      color = "Chromosome"
    ) +
    ggplot2::theme_bw()

}
