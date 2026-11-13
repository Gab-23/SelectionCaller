plot_GC_bias <- function(topTags, gc, n_bins = 20, title_suffix = ""){
  # Diagnostic of GC-content bias in differential accessibility (EDASeq-style).
  # topTags: data.frame with binID, logFC, logCPM (output of process_topTags)
  # gc:      data.frame with binID, GC_content (output of process_GC_content)
  # Returns a ggplot with two facets (logFC, logCPM) of violins by GC-quantile bin.
  # Unbiased data: logFC medians flat around 0 and similar spread across GC bins.

  required_tt <- c("binID", "logFC", "logCPM")
  required_gc <- c("binID", "GC_content")
  if (!all(required_tt %in% colnames(topTags))) stop("plot_GC_bias: topTags must contain ", paste(required_tt, collapse = ", "))
  if (!all(required_gc %in% colnames(gc))) stop("plot_GC_bias: gc must contain ", paste(required_gc, collapse = ", "))

  df <- merge(topTags[, required_tt], gc[, required_gc], by = "binID")              # explicit ID matching, not positional
  df <- df[is.finite(df$logFC) & is.finite(df$logCPM) & is.finite(df$GC_content), ]
  if (nrow(df) < n_bins) stop("plot_GC_bias: fewer matched bins (", nrow(df), ") than GC groups (", n_bins, ")")

  breaks <- unique(stats::quantile(df$GC_content, probs = seq(0, 1, length.out = n_bins + 1)))  # equal-count GC groups
  df$gc_group <- cut(df$GC_content, breaks = breaks, include.lowest = TRUE, dig.lab = 3)

  rho_signed <- stats::cor(df$logFC, df$GC_content, method = "spearman")             # bias in the mean
  rho_abs    <- stats::cor(abs(df$logFC), df$GC_content, method = "spearman")        # bias in the spread
  rho_cpm    <- stats::cor(df$logCPM, df$GC_content, method = "spearman")            # GC effect on depth

  long <- rbind(
    data.frame(gc_group = df$gc_group, value = df$logFC,  metric = "logFC"),
    data.frame(gc_group = df$gc_group, value = df$logCPM, metric = "logCPM")
  )
  long$metric <- factor(long$metric, levels = c("logFC", "logCPM"))
  zero_line <- data.frame(metric = factor("logFC", levels = levels(long$metric)), y = 0)

  p <- ggplot2::ggplot(long, ggplot2::aes(x = gc_group, y = value, color = gc_group)) +
    ggplot2::geom_violin(width = 0.95, scale = "width") +
    ggplot2::geom_boxplot(width = 0.15, color = "grey20", outlier.shape = NA) +
    ggplot2::geom_hline(data = zero_line, ggplot2::aes(yintercept = y), color = "black", lty = 2) +
    ggplot2::facet_wrap(~ metric, ncol = 1, scales = "free_y") +
    ggplot2::scale_color_viridis_d(option = "C", end = 0.9) +
    ggplot2::labs(
      title = paste0("logFC and logCPM by GC-content bin", title_suffix),
      subtitle = sprintf("Spearman with GC: logFC = %.2f (mean bias) | |logFC| = %.2f (spread bias) | logCPM = %.2f | n = %d bins",
                         rho_signed, rho_abs, rho_cpm, nrow(df)),
      x = "GC-content bin (equal-count quantiles)", y = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
                   legend.position = "none",
                   axis.title = ggplot2::element_text(size = 14))

  return(p)
}