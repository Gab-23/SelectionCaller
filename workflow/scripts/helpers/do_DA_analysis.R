do_DA_analysis <- function(counts, library_sizes, mm, gc = NULL, gc_correction = "none", gc_span = 0.3, gc_n_fit = 20000){
# gc_correction: "none" (previous behaviour, unchanged) | "loess" (per-sample GC offset, see compute_GC_offset.R)
# gc:            data.frame with binID, GC_content (required when gc_correction != "none")

if (!gc_correction %in% c("none", "loess")) stop("do_DA_analysis: gc_correction must be 'none' or 'loess', got '", gc_correction, "'")
if (gc_correction != "none" && is.null(gc)) stop("do_DA_analysis: gc_correction = '", gc_correction, "' requires gc")

rownames(counts) <- counts$binID
counts$binID <- NULL

#rownames(counts_sva) <- counts_sva$binID
#counts_sva$binID <- NULL

libsize <- library_sizes$library_size
names(libsize) <- library_sizes$SRA_id

stopifnot(                                                                      # fix: was a bare `==` comparison with the result never checked - both
    "library_sizes and counts are not aligned by sample" =                      # ultimately derive their sample order from list.files() on the same
    all(names(libsize) == colnames(counts))                                    # coverage directory, so this likely already held in practice, but
)                                                                               # nothing was actually enforcing it (same class of issue as the
                                                                                 # binID-alignment check in do_get_counts_analysis.R).

gc_offset <- NULL                                                               # NEW: optional per-sample GC offset
if (gc_correction == "loess") {
    gc_offset <- compute_GC_offset(counts, libsize, gc, span = gc_span, n_fit = gc_n_fit)
}

topTags <- edg2ATAC(counts, mm, libsize, "conditiontumor", gc_offset = gc_offset)
topTags$binID <- rownames(topTags)
rownames(topTags) <- NULL

#manual_topTags <- compute_FC_manually(counts_sva, metadata)
#colnames(manual_topTags) <- c("binID", "NAT", "tumor", "FC", "logFC")
#manual_topTags$PValue <- rep(1, nrow(manual_topTags))
#manual_topTags$FDR <- rep(1, nrow(manual_topTags))

topTags <- process_topTags(topTags)
#manual_topTags <- process_topTags(manual_topTags)
  
p_diagnostics <- plot_diagnostics(topTags)
p_volcano <- plot_volcano(topTags)
p_landscape <- plot_landscape_FC(topTags)
#p_landscape_manual <- plot_landscape_FC(manual_topTags)

#p_comparison <- plot_comparison_FC(topTags, manual_topTags)
p_CNA_comparison <- compare_acc_CNAs(topTags, ampl, del)

out.list <- list(topTags = topTags,
                 p_diagnostics = p_diagnostics,
                 p_volcano = p_volcano,
                 p_landscape = p_landscape,
                 p_CNA_comparison = p_CNA_comparison,
                 gc_offset = gc_offset)                                         # NEW: NULL when gc_correction = "none"

return(out.list)

}