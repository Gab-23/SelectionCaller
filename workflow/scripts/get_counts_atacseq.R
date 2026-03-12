source("workflow/scripts/helpers/process_acc_vals.R")

options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_dir                             <- args$input_dir
output_path_counts                    <- args$output_path_counts
output_path_PCA                       <- args$output_path_PCA

files <- list.files(input_dir, recursive = T, full.names = T)
sample_names <- unlist(lapply(X = files, FUN = function(x){strsplit(x, split = "/")[[1]][[2]]}))

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

count.f.tmm <- process_acc_vals(counts = counts, 
                                minCount_acc = minCounts_atac, 
                                minSamples_acc = minSamples_atac)

pca <- prcomp(t(count.f.tmm), scale. = TRUE)

pca_df <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[,1],
  PC2 = pca$x[,2]
)

pca_plot <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2, label = Sample)) +
                ggplot2::geom_point(size = 3) +
                ggplot2::geom_text(vjust = -0.5) +
                ggplot2::theme_bw() +
                ggplot2::labs(title = "PCA of ATAC-seq accessibility",
                    x = paste0("PC1 (", round(summary(pca)$importance[2,1] * 100, 2), "%)"),
                    y = paste0("PC2 (", round(summary(pca)$importance[2,2] * 100, 2), "%)"))

write.table(x = count.f.tmm,                                                     # write output
            file = output_path_counts, 
            quote = F, 
            sep = "\t", 
            row.names = F, 
            col.names = T)

ggsave(output_path_PCA, plot = p, width = 6, height = 5)