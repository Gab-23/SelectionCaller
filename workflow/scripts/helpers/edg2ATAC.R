edg2ATAC <- function(e, mm, libsize, coef){
  rownames(mm) <- colnames(e)
  dds <- edgeR::calcNormFactors(DGEList(e, lib.size = libsize))
  dds <- edgeR::estimateGLMRobustDisp(dds, mm)
  dds <- edgeR::glmFit(dds, mm)
  res <- as.data.frame(topTags(glmLRT(dds, coef), nrow(e)))
  res
}
