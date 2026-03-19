source("workflow/scripts/helpers/process_acc_vals.R")

options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_dir                             <- args$input_dir
input_path_metadata                   <- args$input_path_metadata
output_path_counts                    <- args$output_path_counts
output_path_PCA                       <- args$output_path_PCA
minSamples_atacseq                    <- args$minSamples_atacseq
minCounts_atacseq                     <- args$minCounts_atacseq

files <- list.files(input_dir, recursive = T, full.names = T)
metadata <- read.table(input_path_metadata, header = T, sep = "\t")

sample_names <- unlist(lapply(X = files, FUN = function(x){strsplit(x, split = "/")[[1]][[5]]}))

cov_list <- lapply(X = seq(length(files)), 
                   FUN = function(idx){
                    file <- files[idx]
                    sample_name <- sample_names[idx]
                    tab <- read.table(file, header = F, sep = "\t")
                    colnames(tab) <- c("chr", "start", "end", "binID", sample_name)
                    to.return <- tab[,c("binID", sample_name)]
                    return(to.return)
                    })

counts <- Reduce(f = function(x,y){merge(x,y)}, 
                 x = cov_list)

library_sizes <- colSums(counts[,-1])

count.f.tmm <- process_acc_vals(counts = counts,
                                minSamples_acc = minSamples_atacseq,
                                minCount_acc = minCounts_atacseq
                                )

pca <- prcomp(t(count.f.tmm), scale. = TRUE)

pca_df <- data.frame(
  SRA_id = rownames(pca$x),
  PC1 = pca$x[,1],
  PC2 = pca$x[,2]
)

pca_df <- merge(pca_df, metadata, by = "SRA_id")

pca_plot <- ggplot2::ggplot(
  pca_df,
  ggplot2::aes(x = PC1, y = PC2, label = SRA_id, color = condition)
) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_text(vjust = -0.5) +
  ggplot2::theme_bw() +
  ggplot2::labs(
    title = "PCA of ATAC-seq accessibility",
    x = paste0("PC1 (", round(summary(pca)$importance[2,1] * 100, 2), "%)"),
    y = paste0("PC2 (", round(summary(pca)$importance[2,2] * 100, 2), "%)")
  )

write.table(x = count.f.tmm,                                                     # write output
            file = output_path_counts, 
            quote = F, 
            sep = "\t", 
            row.names = F, 
            col.names = T)

ggplot2::ggsave(output_path_PCA, plot = pca_plot, width = 6, height = 5)