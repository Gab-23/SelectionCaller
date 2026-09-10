do_DA_analysis <- function(counts, library_sizes, mm){

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

topTags <- edg2ATAC(counts, mm, libsize, "conditiontumor")
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
                 p_CNA_comparison = p_CNA_comparison)

return(out.list)

}