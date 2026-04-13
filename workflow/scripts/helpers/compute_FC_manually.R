compute_FC_manually <- function(counts_sva, metadata){
	
  t_counts_sva <- as.data.frame(t(counts_sva))
  
  bins <- colnames(t_counts_sva)
  samples <- rownames(t_counts_sva)
  
  conditions <- metadata[match(samples, metadata$SRA_id), "condition"]
  
  t_counts_sva$condition <- conditions
  
  manual_FC <- t_counts_sva %>%
                dplyr::group_by(condition) %>%
                dplyr::summarise(dplyr::across(dplyr::where(is.numeric), \(x) median(x, na.rm = TRUE))) %>%
                tidyr::pivot_longer(-condition, names_to = "feature", values_to = "median") %>%
                tidyr::pivot_wider(names_from = condition, values_from = median) %>%
                dplyr::mutate(FC = tumor / NAT, logFC = log2(tumor / NAT))

}
