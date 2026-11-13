compute_GC_offset <- function(counts, libsize, gc, span = 0.3, n_fit = 20000, seed = 1){
  # Per-sample GC-content offset for edgeR (EDASeq/cqn-style within-sample GC correction).
  #
  # counts:  numeric matrix or data.frame, bins x samples, rownames = binID (no binID column)
  # libsize: named numeric vector of library sizes (names = sample IDs, same order as columns)
  # gc:      data.frame with binID, GC_content
  # span:    loess span (fraction of points in each local fit; larger = smoother curve)
  # n_fit:   number of bins used to FIT each loess (random subset, for speed); the curve is
  #          then evaluated on all bins
  #
  # Steps:
  #   1. effective library size = libsize * TMM factor (same as edg2ATAC)
  #   2. per sample: y_ij = log((count_ij + 0.5) / efflib_j)          (natural log, edgeR scale)
  #   3. per sample: robust loess y ~ GC  ->  g_j(GC_i)                (that sample's GC curve)
  #   4. common curve c_i = mean_j g_j(GC_i)                           (kept: shared by all samples)
  #   5. offset_ij = g_j(GC_i) - c_i                                   (removed: sample-specific part)
  # Bins without GC get offset 0. Returns a bins x samples matrix with the dimnames of counts.

  mat <- as.matrix(counts)
  if (is.null(rownames(mat))) stop("compute_GC_offset: counts must have binIDs as rownames")
  if (!all(colnames(mat) == names(libsize))) stop("compute_GC_offset: libsize names must match counts columns, in order")
  if (!all(c("binID", "GC_content") %in% colnames(gc))) stop("compute_GC_offset: gc must contain binID, GC_content")
  if (!(is.numeric(span) && span > 0 && span <= 1)) stop("compute_GC_offset: span must be in (0, 1]")

  gc_vec <- gc$GC_content[match(rownames(mat), gc$binID)]                          # explicit ID matching
  ok <- is.finite(gc_vec)
  if (sum(ok) < 100) stop("compute_GC_offset: fewer than 100 bins with GC content")

  nf <- edgeR::calcNormFactors(edgeR::DGEList(mat, lib.size = libsize))$samples$norm.factors
  efflib <- libsize * nf
  logy <- log(sweep(mat + 0.5, 2, efflib, "/"))

  set.seed(seed)
  fit_idx <- which(ok); if (length(fit_idx) > n_fit) fit_idx <- sort(sample(fit_idx, n_fit))
  gc_eval <- pmin(pmax(gc_vec[ok], min(gc_vec[fit_idx])), max(gc_vec[fit_idx]))   # no extrapolation outside fitted range

  curves <- matrix(0, nrow(mat), ncol(mat), dimnames = dimnames(mat))
  for (j in seq_len(ncol(mat))) {
    d <- data.frame(y = logy[fit_idx, j], x = gc_vec[fit_idx])
    lo <- stats::loess(y ~ x, data = d, span = span, degree = 1, family = "symmetric")  # robust to DA bins/outliers
    curves[ok, j] <- stats::predict(lo, newdata = data.frame(x = gc_eval))
  }
  offset <- curves - rowMeans(curves)
  offset[!ok, ] <- 0
  offset
}