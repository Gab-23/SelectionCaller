options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_DARs                               <- args$input_path_DARs
input_path_centromeres                        <- args$input_path_centromeres
input_path_telomeres                          <- args$input_path_telomeres
input_path_coding_mutations                   <- args$input_path_coding_mutations
input_path_noncoding_mutations                <- args$input_path_noncoding_mutations
input_path_dist_to_essential                  <- args$input_path_dist_to_essential
input_path_epigenomes                         <- args$input_path_epigenomes
input_path_num_fragile_sites                  <- args$input_path_num_fragile_sites
input_path_dist_fragile_sites                 <- args$input_path_dist_fragile_sites
input_path_GC_content                         <- args$input_path_GC_content
input_path_ampl_frequency                     <- args$input_path_ampl_frequency
input_path_del_frequency                      <- args$input_path_del_frequency
input_path_num_overlaps_protein_coding        <- args$input_path_num_overlaps_protein_coding
input_path_num_overlaps_lncRNA                <- args$input_path_num_overlaps_lncRNA
input_path_num_overlaps_miRNA                 <- args$input_path_num_overlaps_miRNA
input_path_num_overlaps_rRNA                  <- args$input_path_num_overlaps_rRNA
input_path_dist_protein_coding                <- args$input_path_dist_protein_coding
input_path_dist_lncRNA                        <- args$input_path_dist_lncRNA
input_path_dist_miRNA                         <- args$input_path_dist_miRNA
input_path_dist_rRNA                          <- args$input_path_dist_rRNA
input_path_TSG                                <- args$input_path_TSG
input_path_OG                                 <- args$input_path_OG
input_path_repliseq                           <- args$input_path_repliseq
output_path_dataset                           <- args$output_path_dataset
output_path_corr_plot                         <- args$output_path_corr
output_path_NA_plot                           <- args$output_path_NA_plot
output_path_boxplots                          <- args$output_path_boxplots
#output_path_PCA_plot                         <- args$output_path_PCA_plot

`%>%` <- magrittr::`%>%`

input_paths <- lapply(X = grep(pattern = "input_path", 
                               x = ls(), 
                               value = TRUE), 
                      FUN = get)

data.frames <- lapply(X = input_paths, 
                      FUN = function(input_path){
                        to.return <- read.delim(file = input_path,
                                                header = T, 
                                                sep = "\t",
                                                check.names = F
                                                )
                        return(to.return)
                      })

dataset <- Reduce(x = data.frames, 
                  f = function(x,y){
                         merge(x, y, 
                               by="binID", 
                               sort = F, 
                               all = T)
                        }
                    )

idx_cols_to_impute <- grep(pattern = "length_(non)?coding", x = colnames(dataset))
dataset[idx_cols_to_impute][is.na(dataset[idx_cols_to_impute])] <- 0

na_fraction <- data.frame(
  variable = names(dataset[,-1]),
  na_frac = sapply(dataset[,-1], function(x) mean(is.na(x)))
)

long_data <- tidyr::pivot_longer(
  dataset,
  cols = -binID,
  names_to = "variable",
  values_to = "value"
)

cor.mat <- cor(dataset[,-1], use = "pairwise.complete.obs", method = "pearson")
cor.df <- as.data.frame(as.table(cor.mat))
names(cor.df) <- c("Var1", "Var2", "Correlation")

#pca <- prcomp(dataset[,-1], scale. = TRUE)
#pca_variance_vec <- summary(pca)$importance[2,]
#pca_df <- data.frame(SRA_id = rownames(pca$x), PC1 = pca$x[,1], PC2 = pca$x[,2])

NA_frac_plot <- ggplot2::ggplot(
  na_fraction,
  ggplot2::aes(x = stats::reorder(variable, -na_frac), y = na_frac)
) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::theme_bw() +
  ggplot2::labs(x = "Variable", y = "Fraction of NA") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))

boxplots <- ggplot2::ggplot(long_data, ggplot2::aes(x = variable, y = value)) +
  ggplot2::geom_boxplot() +
  ggplot2::theme_bw() +
  ggplot2::labs(x = "Variable", y = "Value") +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))

heatmap <- ggplot2::ggplot(cor.df, ggplot2::aes(x = Var1, y = Var2, fill = Correlation)) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    axis.title = ggplot2::element_blank()
  )

write.table(x = dataset, file = output_path, sep = "\t", quote = F, col.names = T, row.names = F)

ggplot2::ggsave(output_path_NA_plot, NA_frac_plot, width = 12, height = 7)
ggplot2::ggsave(output_path_boxplots, boxplots, width = 12, height = 7)
ggplot2::ggsave(output_path_corr, heatmap, width = 12, height = 12)