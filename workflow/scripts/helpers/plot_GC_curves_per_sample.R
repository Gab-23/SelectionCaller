plot_GC_curves_per_sample <- function(counts, library_sizes, metadata, gc, n_bins = 20, title_suffix = "", gc_offset = NULL){
  # Per-sample GC-content curves (EDASeq/cqn-style diagnostic).
  # counts:        data.frame with a binID column + one column per sample (SRA_id), as read in get_DARs.R
  # library_sizes: data.frame with SRA_id, library_size
  # metadata:      data.frame with SRA_id, condition, study_id
  # gc:            data.frame with binID, GC_content
  # gc_offset:     optional bins x samples matrix from compute_GC_offset() (natural log);
  #                if given, curves are drawn AFTER subtracting it (to check the correction)
  # Top row: median logCPM per GC bin, one line per sample.
  # Bottom row: the same curve minus the across-sample mean curve (what a per-sample GC
  # correction would remove). Columns: study_id; colour: condition.
  # Subtitle: R2 of the per-sample GC slope explained by study vs by condition
  # (slopes grouping by study -> technical; by condition within study -> possibly biological).

  if (!"binID" %in% colnames(counts)) stop("plot_GC_curves_per_sample: counts must contain a binID column")
  if (!all(c("SRA_id", "library_size") %in% colnames(library_sizes))) stop("plot_GC_curves_per_sample: library_sizes must contain SRA_id, library_size")
  if (!all(c("SRA_id", "condition", "study_id") %in% colnames(metadata))) stop("plot_GC_curves_per_sample: metadata must contain SRA_id, condition, study_id")
  if (!all(c("binID", "GC_content") %in% colnames(gc))) stop("plot_GC_curves_per_sample: gc must contain binID, GC_content")

  mat <- as.matrix(counts[, setdiff(colnames(counts), "binID"), drop = FALSE])
  rownames(mat) <- counts$binID
  samples <- colnames(mat)
  if (!all(samples %in% library_sizes$SRA_id)) stop("plot_GC_curves_per_sample: samples without library size: ",
                                                    paste(setdiff(samples, library_sizes$SRA_id), collapse = ", "))
  if (!all(samples %in% metadata$SRA_id)) stop("plot_GC_curves_per_sample: samples without metadata: ",
                                               paste(setdiff(samples, metadata$SRA_id), collapse = ", "))

  gc_vec <- gc$GC_content[match(rownames(mat), gc$binID)]                         # explicit ID matching
  keep <- is.finite(gc_vec)
  mat <- mat[keep, , drop = FALSE]; gc_vec <- gc_vec[keep]
  if (nrow(mat) < n_bins) stop("plot_GC_curves_per_sample: fewer bins with GC (", nrow(mat), ") than GC groups (", n_bins, ")")

  libsize <- library_sizes$library_size[match(samples, library_sizes$SRA_id)]
  dge <- edgeR::calcNormFactors(edgeR::DGEList(mat, lib.size = libsize))            # same normalisation as edg2ATAC
  logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 2)
  if (!is.null(gc_offset)) {                                                         # corrected curves
    if (!all(rownames(mat) %in% rownames(gc_offset)) || !all(samples %in% colnames(gc_offset)))
      stop("plot_GC_curves_per_sample: gc_offset does not cover all bins/samples")
    logcpm <- logcpm - gc_offset[rownames(mat), samples] / log(2)                    # natural log -> log2
  }

  breaks <- unique(stats::quantile(gc_vec, probs = seq(0, 1, length.out = n_bins + 1)))
  grp <- cut(gc_vec, breaks = breaks, include.lowest = TRUE)
  gc_mid <- tapply(gc_vec, grp, stats::median)                                       # x position of each GC group
  curves <- apply(logcpm, 2, function(v) tapply(v, grp, stats::median))              # GC groups x samples
  deviation <- curves - rowMeans(curves)                                             # sample-specific part

  meta <- metadata[match(samples, metadata$SRA_id), c("SRA_id", "condition", "study_id")]
  to_long <- function(m, panel) data.frame(
    SRA_id = rep(colnames(m), each = nrow(m)), GC = rep(as.numeric(gc_mid), times = ncol(m)),
    value = as.vector(m), panel = panel)
  long <- rbind(to_long(curves, "median logCPM"), to_long(deviation, "deviation from mean curve"))
  long <- merge(long, meta, by = "SRA_id")
  long$panel <- factor(long$panel, levels = c("median logCPM", "deviation from mean curve"))

  # per-sample slope of the deviation curve, and how much of it study vs condition explain
  slopes <- data.frame(SRA_id = colnames(deviation),
                       slope = apply(deviation, 2, function(d) stats::coef(stats::lm(d ~ as.numeric(gc_mid)))[2]))
  slopes <- merge(slopes, meta, by = "SRA_id")
  r2 <- function(f) if (length(unique(slopes[[all.vars(f)[2]]])) > 1) summary(stats::lm(f, data = slopes))$r.squared else NA_real_
  r2_study <- r2(slope ~ study_id)
  r2_cond  <- r2(slope ~ condition)
  r2_both  <- if (length(unique(slopes$study_id)) > 1 && length(unique(slopes$condition)) > 1)
                summary(stats::lm(slope ~ study_id + condition, data = slopes))$r.squared else NA_real_

  p <- ggplot2::ggplot(long, ggplot2::aes(x = GC, y = value, group = SRA_id, color = condition)) +
    ggplot2::geom_line(alpha = 0.7) +
    ggplot2::facet_grid(panel ~ study_id, scales = "free_y") +
    ggplot2::labs(
      title = paste0("Per-sample GC curves", title_suffix),
      subtitle = sprintf("R2 of per-sample GC slope: ~ study %.2f | ~ condition %.2f | ~ study + condition %.2f | %d samples, %d bins",
                         r2_study, r2_cond, r2_both, ncol(mat), nrow(mat)),
      x = "GC content (median of GC-quantile group)", y = NULL, color = "condition") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title = ggplot2::element_text(size = 14), legend.position = "bottom")

  attr(p, "slopes") <- slopes                                                        # for tests / downstream inspection
  return(p)
}