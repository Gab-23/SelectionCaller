options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_topTags_1                     <- args$input_path_topTags_1 #"/hpcnfs/scratch/MS/ieo7429/SelectionCaller/DARs_n, header = T, sep = "\t"
input_path_topTags_2	                 <- args$input_path_topTags_2 #"/hpcnfs/scratch/MS/ieo7429/SelectionCaller/optional_DARs, header = T, sep = "\t"
suffix_1                                 <- args$suffix_1
suffix_2                                 <- args$suffix_2

output_path_confusion_plot_significant   <- args$output_path_DARs
output_path_confusion_plot_sign		     <- args$output_path_manual_DARs
output_path_scatter_plot	             <- args$output_path_diagnostics_plot

topTags_1 <- read.table(input_path_topTags_1, header = T, sep = "\t")
topTags_2 <- read.table(input_path_topTags_2, header = T, sep = "\t")

topTags_1 <- topTags_1[, c("binID", "logFC", "logCPM", "PValue", "FDR", "significant")]
topTags_2 <- topTags_2[, c("binID", "logFC", "logCPM", "PValue", "FDR", "significant")]

topTags_1$sign <- sign(topTags_1$logFC)
topTags_2$sign <- sign(topTags_2$logFC)

merged_topTags <- merge(x = topTags_1, y = topTags_2, by = "binID", suffixes = c(suffix_1, suffix_2))

conf_mat_significant <- table(
merged_topTags[[paste0("significant", suffix_1)]],
merged_topTags[[paste0("significant", suffix_2)]]
)

conf_df_significant <- as.data.frame(conf_mat_significant)
colnames(conf_df_significant) <- c("Set1", "Set2", "Freq")

confusion_plot_significant <- ggplot2::ggplot(
								conf_df_significant,
								ggplot2::aes(x = Set2, y = Set1, fill = Freq)
								) +
								ggplot2::geom_tile(color = "black") +
								ggplot2::geom_text(ggplot2::aes(label = Freq)) +
								ggplot2::scale_fill_gradient(low = "white", high = "steelblue") +
								ggplot2::coord_fixed() +
								ggplot2::labs(
								  title = "Confusion Matrix of Significant DARs",
								  x = paste0("significant", suffix_2),
								  y = paste0("significant", suffix_1)
								) +
								ggplot2::theme_bw()

conf_mat_sign <- table(
merged_topTags[[paste0("sign", suffix_1)]],
merged_topTags[[paste0("sign", suffix_2)]]
)

conf_df_sign <- as.data.frame(conf_mat_sign)
colnames(conf_df_sign) <- c("Set1", "Set2", "Freq")

confusion_plot_sign <- ggplot2::ggplot(
							conf_df_sign,
							ggplot2::aes(x = Set2, y = Set1, fill = Freq)
							) +
							ggplot2::geom_tile(color = "black") +
							ggplot2::geom_text(ggplot2::aes(label = Freq)) +
							ggplot2::scale_fill_gradient(low = "white", high = "darkred") +
							ggplot2::coord_fixed() +
							ggplot2::labs(
							  title = "Confusion Matrix of logFC Sign",
							  x = paste0("sign", suffix_2),
							  y = paste0("sign", suffix_1)
							) +
							ggplot2::theme_bw()

merged_topTags$agreement <- merged_topTags[[paste0("significant", suffix_1)]] ==
                          merged_topTags[[paste0("significant", suffix_2)]]

scatter_plot <- ggplot2::ggplot(
					merged_topTags,
					ggplot2::aes(
					  x = .data[[paste0("logFC", suffix_1)]],
					  y = .data[[paste0("logFC", suffix_2)]],
					  color = agreement
					)
					) +
					ggplot2::geom_point(alpha = 0.6) +
					ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
					ggplot2::labs(
					  title = "logFC Comparison",
					  x = paste0("logFC", suffix_1),
					  y = paste0("logFC", suffix_2)
					) +
					ggplot2::theme_bw()

ggplot2::ggsave(output_path_confusion_plot_significant, confusion_plot_significant, width = 7, height = 7)
ggplot2::ggsave(output_path_confusion_plot_sign, confusion_plot_sign, width = 7, height = 7)
ggplot2::ggsave(output_path_scatter_plot, scatter_plot, width = 7, height = 7)
