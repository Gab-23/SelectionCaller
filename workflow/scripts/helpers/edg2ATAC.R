edg2ATAC <- function(e, mm, libsize, coef, gc_offset = NULL){
  # gc_offset (optional): bins x samples matrix from compute_GC_offset(), natural-log scale,
  # same dimnames as e. If given, the GLM offset becomes log(lib.size * TMM factor) + gc_offset.
  # NOTE: in edgeR an explicit offset REPLACES lib.size/norm.factors in the GLM, which is why
  # the library-size term is added back here. glmFit and all dispersion estimators use it
  # (checked on edgeR 4.0.16); logCPM in topTags does not (it is computed from lib.size only).
  rownames(mm) <- colnames(e)
  dds <- edgeR::calcNormFactors(DGEList(e, lib.size = libsize))
  if (!is.null(gc_offset)) {
    if (!identical(dim(gc_offset), dim(dds$counts)) ||
        !identical(rownames(gc_offset), rownames(dds$counts)) ||
        !identical(colnames(gc_offset), colnames(dds$counts)))
      stop("edg2ATAC: gc_offset must have the same dimensions and dimnames as the counts")
    base <- matrix(log(dds$samples$lib.size * dds$samples$norm.factors),
                   nrow = nrow(dds$counts), ncol = ncol(dds$counts), byrow = TRUE)
    dds$offset <- base + gc_offset
  }
  dds <- edgeR::estimateGLMRobustDisp(dds, mm)
  dds <- edgeR::glmFit(dds, mm)
  res <- as.data.frame(topTags(glmLRT(dds, coef), nrow(e)))
  res
}