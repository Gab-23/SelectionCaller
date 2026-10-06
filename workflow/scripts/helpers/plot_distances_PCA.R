plot_distances_PCA <- function(pca, pca_sva, metadata){

studies <- unique(metadata$study_id)
conditions <- unique(metadata$condition)

distance_matrix <- as.data.frame(as.matrix(dist(pca$x)))
distance_matrix_sva <- as.data.frame(as.matrix(dist(pca_sva$x)))

inter_study_distance <- distance_matrix[rownames(distance_matrix) %in% metadata[metadata$study_id %in% studies[1], "SRA_id"], 
									    colnames(distance_matrix) %in% metadata[metadata$study_id %in% studies[2], "SRA_id"]]
								   
inter_condition_distance <- distance_matrix[rownames(distance_matrix) %in% metadata[metadata$condition %in% conditions[1], "SRA_id"], 
									        colnames(distance_matrix) %in% metadata[metadata$condition %in% conditions[2], "SRA_id"]]
									        
inter_study_distance_sva <- distance_matrix_sva[rownames(distance_matrix_sva) %in% metadata[metadata$study_id %in% studies[1], "SRA_id"], 
									            colnames(distance_matrix_sva) %in% metadata[metadata$study_id %in% studies[2], "SRA_id"]]
								   
inter_condition_distance_sva <- distance_matrix_sva[rownames(distance_matrix_sva) %in% metadata[metadata$condition %in% conditions[1], "SRA_id"], 
									                colnames(distance_matrix_sva) %in% metadata[metadata$condition %in% conditions[2], "SRA_id"]]

inter_study_distance_vec <- unlist(inter_study_distance)		   
inter_condition_distance_vec <- unlist(inter_condition_distance)

inter_study_distance_sva_vec <- unlist(inter_study_distance_sva)		   
inter_condition_distance_sva_vec <- unlist(inter_condition_distance_sva)

# matched donors: NAT vs tumour distance of the SAME donor, paired explicitly by donor_id
# (fix: previous match()/diag() logic misaligned rows whenever a donor had only one condition)
get_matched_distances <- function(inter_condition_distance){
  donor_rows <- metadata[match(rownames(inter_condition_distance), metadata$SRA_id), "donor_id"]
  donor_cols <- metadata[match(colnames(inter_condition_distance), metadata$SRA_id), "donor_id"]
  common_donor_ids <- intersect(donor_rows, donor_cols)
  vapply(common_donor_ids, function(d){
    stopifnot(sum(donor_rows == d) == 1, sum(donor_cols == d) == 1)
    inter_condition_distance[donor_rows == d, donor_cols == d]
  }, numeric(1))
}

matched_distances     <- get_matched_distances(inter_condition_distance)
matched_distances_sva <- get_matched_distances(inter_condition_distance_sva)


distance_df <- data.frame(
distance = c(
  unlist(inter_study_distance), unlist(inter_condition_distance), matched_distances,
  unlist(inter_study_distance_sva), unlist(inter_condition_distance_sva), matched_distances_sva
),
type = factor(c(
  rep("Inter-study", length(unlist(inter_study_distance))),
  rep("Inter-condition", length(unlist(inter_condition_distance))),
  rep("Matched donors", length(matched_distances)),
  rep("Inter-study", length(unlist(inter_study_distance_sva))),
  rep("Inter-condition", length(unlist(inter_condition_distance_sva))),
  rep("Matched donors", length(matched_distances_sva))
),
levels = c("Inter-study", "Inter-condition", "Matched donors")),
method = rep(c("Raw", "SVA"),
             times = c(
               length(unlist(inter_study_distance)) + length(unlist(inter_condition_distance)) + length(matched_distances),
               length(unlist(inter_study_distance_sva)) + length(unlist(inter_condition_distance_sva)) + length(matched_distances_sva)
             ))
)
  

p <- ggplot2::ggplot(distance_df, ggplot2::aes(x = type, y = distance, fill = method)) +
ggplot2::geom_boxplot(position = ggplot2::position_dodge(width = 0.8), alpha = 0.7, outlier.shape = NA) +
ggplot2::theme_bw() +
ggplot2::labs(title = "Distance Distributions Before and After SVA", x = "", y = "Distance") +
ggplot2::scale_fill_manual(values = c("Raw" = "#F8766D", "SVA" = "#00BFC4")) +
ggplot2::theme(legend.title = ggplot2::element_blank())

}