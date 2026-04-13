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

corresponding_donor_id_rows <- metadata[match(rownames(inter_condition_distance), metadata$SRA_id), "donor_id"]
corresponding_donor_id_cols <- metadata[match(colnames(inter_condition_distance), metadata$SRA_id), "donor_id"]
common_donor_ids <- intersect(corresponding_donor_id_rows, corresponding_donor_id_cols)

matching_donor_ids <- match(corresponding_donor_id_rows[corresponding_donor_id_rows %in% common_donor_ids], 
			    corresponding_donor_id_cols[corresponding_donor_id_cols %in% common_donor_ids])

corresponding_donor_id_rows_sva <- metadata[match(rownames(inter_condition_distance_sva), metadata$SRA_id), "donor_id"]
corresponding_donor_id_cols_sva <- metadata[match(colnames(inter_condition_distance_sva), metadata$SRA_id), "donor_id"]
common_donor_ids_sva <- intersect(corresponding_donor_id_rows_sva, corresponding_donor_id_cols_sva)

matching_donor_ids_sva <- match(corresponding_donor_id_rows_sva[corresponding_donor_id_rows_sva %in% common_donor_ids_sva],
                                corresponding_donor_id_cols_sva[corresponding_donor_id_cols_sva %in% common_donor_ids_sva])

inter_condition_distance_matched <- inter_condition_distance[, matching_donor_ids]
matched_distances <- diag(as.matrix(inter_condition_distance_matched))

inter_condition_distance_matched_sva <- inter_condition_distance_sva[, matching_donor_ids_sva]
matched_distances_sva <- diag(as.matrix(inter_condition_distance_matched_sva))


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
