source("workflow/scripts/helpers/find_neighbors.R")
source("workflow/scripts/helpers/plot_moran_autocorrelation.R")

library(spdep)

options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

plots_dir		              <- args$plots_dir
input_path_dataset                    <- args$input_path_dataset
input_path_backbone                   <- args$input_path_backbone
K				      <- as.integer(args$K)

output_path_univariate_table	      <- args$output_path_univariate_table
output_path_plot_univariate           <- args$output_path_plot_univariate

dataset <- read.delim(input_path_dataset, check.names = F, header = T, sep = "\t")
colnames(dataset) <- gsub(pattern = "/", replacement = "_", x = colnames(dataset))

backbone <- read.table(input_path_backbone, header = F, sep = "\t")
colnames(backbone) <- c("chr", "start", "end", "binID", "annot")

nb <- lapply(X = 1:nrow(backbone), 
	   FUN = function(idx){
	           neighbors <- find_neighbors(backbone, idx, K = K)
	   }
	  )

#not_null_mask <- which(!(unlist(lapply(nb, is.null))))
#backbone.f <- backbone[not_null_mask,]
#dataset.f <- dataset[dataset$binID %in% backbone.f$binID, ]
#nb.f <- nb[not_null_mask]

class(nb) <- "nb"
listw <- nb2listw(nb)

moran_results <- do.call(rbind,
                         lapply(X = (1:ncol(dataset))[-1], 
                                FUN = function(idx){
				      feature_name <- colnames(dataset)[idx]
                                      moran_result <- moran.test(log1p(dataset[,idx]), listw = listw, randomisation = T)
                                      pval <- round(moran_result$p.value, 10)
			              statistic <- moran_result$statistic
			              df.to.return <- data.frame(statistic = statistic, pval = pval, row.names = feature_name)
                                }))

moran_results$FDR <- p.adjust(moran_results$pval, "fdr")
moran_results$state <- rownames(moran_results)
rownames(moran_results) <- NULL

moran_results$state <- factor(
  moran_results$state,
  levels = moran_results$state[order(moran_results$statistic, decreasing = TRUE)]
)

p <- ggplot2::ggplot(moran_results, ggplot2::aes(x = state, y = statistic)) +
     ggplot2::geom_bar(stat = "identity") +
     ggplot2::coord_flip() +
     ggplot2::labs(x = "State", y = "Statistic", title = "Moran Results Statistic") +
     ggplot2::theme_bw()

lapply(X = (1:ncol(dataset))[-1],
       FUN = function(idx){
    	  
    	  col <- dataset[,idx]
    	  col_name <- colnames(dataset)[idx]
    	  
          plot_data <- moran.plot(x = log1p(col), listw = listw, labels = dataset$binID)
          plot_data.s <- plot_data[,c("x", "wx", "labels")]
	  colnames(plot_data.s) <- c("x", "y", "binID")
	  
 	  plot_data.s <- merge(x = plot_data.s, y = backbone[,c("binID", "annot")], by = "binID")
 
	  plot_autocorrelation <- plot_moran_autocorrelation(plot_data.s, K = K)
          
	  ggplot2::ggsave(paste0(plots_dir, K, "_autocorrelation_annotated_", col_name, ".pdf"), plot_autocorrelation, width = 7, height = 7)

    }
  )

ggplot2::ggsave(output_path_plot_univariate, p, width = 7, height = 7)
write.table(x = moran_results, file = output_path_univariate_table, sep = "\t", col.names = T, row.names = F, quote = F)
