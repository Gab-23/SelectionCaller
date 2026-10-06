# Panel: PCA before/after SVA (donors joined) + per-pair distance change in full PCA space.
# Returns ONE object (gtable): save it like any other plot with ggplot2::ggsave(path, plot = ..., width = 13, height = 10).

plot_PCA_SVA_distance_panel <- function(pca, pca_sva, metadata){

  pcas <- list("Before SVA" = pca, "After SVA" = pca_sva)
  samples <- rownames(pca$x)
  stopifnot(identical(samples, rownames(pca_sva$x)), all(samples %in% metadata$SRA_id))
  m <- metadata[match(samples, metadata$SRA_id), ]

  # ---- PCA coordinates (PC1/PC2) ----
  pca_df <- do.call(rbind, lapply(names(pcas), function(k){
    v <- round(summary(pcas[[k]])$importance[2, 1:2] * 100, 1)
    data.frame(SRA_id = samples, PC1 = pcas[[k]]$x[, 1], PC2 = pcas[[k]]$x[, 2],
               method = paste0(k, "  (PC1 ", v[1], "%, PC2 ", v[2], "%)"))
  }))
  pca_df <- merge(pca_df, m[, c("SRA_id", "donor_id", "condition", "study_id")], by = "SRA_id")
  pca_df$method    <- factor(pca_df$method, levels = unique(pca_df$method[order(!grepl("^Before", pca_df$method))]))
  pca_df$condition <- factor(pca_df$condition, levels = c("NAT", "tumor"), labels = c("NAT", "Tumour"))

  # segments joining NAT and tumour of the same donor
  seg_df <- merge(pca_df[pca_df$condition == "NAT",    c("donor_id", "method", "PC1", "PC2")],
                  pca_df[pca_df$condition == "Tumour", c("donor_id", "method", "PC1", "PC2")],
                  by = c("donor_id", "method"), suffixes = c("", "_end"))

  # ---- pairwise distances in full PCA space (all PCs) ----
  pair_idx <- t(utils::combn(length(samples), 2))
  dists <- lapply(pcas, function(p) as.matrix(stats::dist(p$x))[samples, samples][pair_idx])
  same_donor <- m$donor_id[pair_idx[, 1]]  == m$donor_id[pair_idx[, 2]]
  same_study <- m$study_id[pair_idx[, 1]]  == m$study_id[pair_idx[, 2]]
  same_cond  <- m$condition[pair_idx[, 1]] == m$condition[pair_idx[, 2]]
  cat_levels <- c("Same donor\n(NAT vs tumour)", "Same study,\nsame condition",
                  "Different study,\nsame condition", "Same study,\ndifferent condition")
  category <- ifelse(same_donor, cat_levels[1],
              ifelse(same_cond &  same_study, cat_levels[2],
              ifelse(same_cond & !same_study, cat_levels[3],
              ifelse(!same_cond & same_study, cat_levels[4], NA))))   # diff study AND condition: confounded -> dropped

  dist_df <- data.frame(pair = seq_along(category), category = category,
                        Before = dists[["Before SVA"]], After = dists[["After SVA"]])
  dist_df <- dist_df[!is.na(dist_df$category), ]
  dist_df$category   <- factor(dist_df$category, levels = cat_levels)
  dist_df$log2_ratio <- log2(dist_df$After / dist_df$Before)

  dist_long <- rbind(data.frame(dist_df[, c("pair", "category")], method = "Before SVA", distance = dist_df$Before),
                     data.frame(dist_df[, c("pair", "category")], method = "After SVA",  distance = dist_df$After))
  dist_long$method <- factor(dist_long$method, levels = c("Before SVA", "After SVA"))

  # ---- A: PCA before / after, donors joined ----
  p_pca <- ggplot2::ggplot(pca_df, ggplot2::aes(PC1, PC2)) +
    ggplot2::geom_segment(data = seg_df, ggplot2::aes(xend = PC1_end, yend = PC2_end), color = "grey60") +
    ggplot2::geom_point(ggplot2::aes(color = study_id, shape = condition), size = 3) +
    ggplot2::facet_wrap(~method, scales = "free") +
    ggplot2::labs(color = "Study", shape = "Condition", tag = "A") +
    ggplot2::theme_bw()

  # ---- B: per-pair distance before -> after ----
  p_dist <- ggplot2::ggplot(dist_long, ggplot2::aes(x = method, y = distance)) +
    ggplot2::geom_line(ggplot2::aes(group = pair), color = "grey70", alpha = 0.4) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = method), alpha = 0.8, outlier.shape = NA, width = 0.5) +
    ggplot2::scale_fill_manual(values = c("Before SVA" = "#F8766D", "After SVA" = "#00BFC4")) +
    ggplot2::facet_wrap(~category, nrow = 1) +
    ggplot2::labs(x = "", y = "Distance in PCA space", tag = "B") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")

  # ---- C: relative change per pair ----
  p_ratio <- ggplot2::ggplot(dist_df, ggplot2::aes(x = category, y = log2_ratio)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_boxplot(fill = "grey85", outlier.shape = NA, width = 0.5) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.15, seed = 1), size = 0.6, alpha = 0.4) +
    ggplot2::labs(x = "", y = "log2(distance after / before SVA)", tag = "C") +
    ggplot2::theme_bw()

  # ---- assemble A (top, full width), B and C (bottom) into one gtable (gtable ships with ggplot2) ----
  panel <- gtable::gtable(widths = grid::unit(c(1.7, 1), "null"), heights = grid::unit(c(1, 1), "null"))
  panel <- gtable::gtable_add_grob(panel, ggplot2::ggplotGrob(p_pca),   t = 1, l = 1, r = 2)
  panel <- gtable::gtable_add_grob(panel, ggplot2::ggplotGrob(p_dist),  t = 2, l = 1)
  panel <- gtable::gtable_add_grob(panel, ggplot2::ggplotGrob(p_ratio), t = 2, l = 2)
  panel
}