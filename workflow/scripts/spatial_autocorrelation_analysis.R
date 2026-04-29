# do this on data before merging

source("workflow/scripts/helpers/find_neighbors.R")
source("workflow/scripts/helpers/plot_moran_autocorrelation.R")
source("workflow/scripts/helpers/visualize_local_moran.R")
source("workflow/scripts/helpers/GoEnrichment.R")
source("workflow/scripts/helpers/plotGO.R")
source("workflow/scripts/helpers/pva.R")

library(spdep)
library(parallel)

options(scipen=999)                                                               # unable scientific notation  

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

go_dir                                <- args$go_dir
plots_dir		                          <- args$plots_dir
input_path_dataset                    <- args$input_path_dataset
input_path_backbone                   <- args$input_path_backbone
K				                              <- as.integer(args$K)
C                                     <- as.integer(args$C)
FDR_thr                               <- 0.05
ID                                    <- "SYMBOL"

invisible(spdep::set.coresOption(C))

output_path_univariate_table	        <- args$output_path_univariate_table
output_path_plot_univariate           <- args$output_path_plot_univariate
output_path_bool_autocorrelation      <- args$output_path_bool_autocorrelation

`%>%` <- magrittr::`%>%`

#input_path_datasets <- list.files(path = input_dir_datasets, full.names = T, recursive = T)
#for(input_path_dataset in input_path_datasets){
#print(input_path_dataset)

dataset <- read.delim(input_path_dataset, check.names = F, header = T, sep = "\t")
colnames(dataset) <- gsub(pattern = "/", replacement = "_", x = colnames(dataset)) # clean colnames of chromatin states. 20_ZNF/Rpts -> 20_ZNF_Rpts

backbone <- read.table(input_path_backbone, header = F, sep = "\t")
colnames(backbone) <- c("chr", "start", "end", "binID")# , "annot")

# do it chromosome specific?
dataset <- dataset[startsWith(prefix = "1_", x = dataset$binID),]
dataset <- dataset[dataset$binID %in% backbone$binID,]
backbone <- backbone[backbone$binID %in% dataset$binID,]

nb <- parallel::mclapply(
     mc.cores = C,
     X = 1:nrow(backbone), 
	   FUN = function(idx){
	           neighbors <- find_neighbors(backbone, idx, K = K, mode = "spatial_analysis")
	   }
	  )

class(nb) <- "nb"
listw <- nb2listw(nb)

#global_moran_results <- do.call(rbind,
#                         lapply(X = (1:ncol(dataset))[-1], 
#                                FUN = function(idx){
#				                              feature_name <- colnames(dataset)[idx]
#                                      moran_result <- moran.test(log1p(dataset[,idx]), listw = listw, randomisation = T)
#                                      pval <- round(moran_result$p.value, 10)
#			                                statistic <- moran_result$statistic
#			                                df.to.return <- data.frame(statistic = statistic, pval = pval, row.names = feature_name)
#                                }))
#
#global_moran_results$FDR <- p.adjust(global_moran_results$pval, "fdr")
#global_moran_results$state <- rownames(global_moran_results)
#rownames(global_moran_results) <- NULL
#
#global_moran_results$state <- factor(
#  global_moran_results$state,
#  levels = global_moran_results$state[order(global_moran_results$statistic, decreasing = TRUE)]
#)
#
#p <- ggplot2::ggplot(global_moran_results, ggplot2::aes(x = state, y = statistic)) +
#     ggplot2::geom_bar(stat = "identity") +
#     ggplot2::coord_flip() +
#     ggplot2::labs(x = "State", y = "Statistic", title = "Moran Results Statistic") +
#     ggplot2::theme_bw()
#
#valid_symbols <- keys(org.Hs.eg.db, keytype = ID)
#
#lapply(X = (1:ncol(dataset))[-1],
#       FUN = function(idx){
#    	  col <- dataset[,idx]
#    	  col_name <- colnames(dataset)[idx]
#
#        plot_data <- moran.plot(x = log1p(col), listw = listw, labels = dataset$binID)
#        plot_data.s <- plot_data[,c("x", "wx", "labels")]
#    	  colnames(plot_data.s) <- c("x", "y", "binID")
#        plot_data.s <- merge(x = plot_data.s, y = backbone[,c("binID", "annot")], by = "binID")
#
#        clus <- NbClust::NbClust(data = plot_data.s[plot_data.s$x > 0, c("y")], 
#                                 distance = "euclidean", 
#                                 method = "kmeans", 
#                                 index = "all",
#                                 min.nc = 3,
#                                 max.nc = 7
#                )
#
#        plot_data.s$cluster <- 999
#        plot_data.s[plot_data.s$x > 0,]$cluster <- clus$Best.partition
#
#        clusters <- unique(plot_data.s$cluster)
#        allGenes <- unlist(lapply(X = 1:nrow(plot_data.s), 
#                                  FUN = function(idx){
#                                    genes <- strsplit(x = plot_data.s[idx, "annot"], split = ",")[[1]]
#                                    }
#                                  )
#                                )
#        allGenes <- allGenes[allGenes %in% valid_symbols]
#        
#        for (cluster in clusters){                                      # For cluster in clusters
#          cluster_data <- plot_data.s[plot_data.s$cluster == cluster, ] # subset dataset
#          clusterGenes <- unlist(lapply(X = 1:nrow(cluster_data), 
#                                  FUN = function(idx){
#                                    genes <- strsplit(x = cluster_data[idx, "annot"], split = ",")[[1]]
#                                    }
#                                  )
#                                )
#          clusterGenes <- clusterGenes[clusterGenes %in% valid_symbols]
#
#          go <- go.enrichment(allGenes = unique(allGenes), 
#                              deGenes = unique(clusterGenes), 
#                              ID = ID, 
#                              cutoff = 1)
#
#          for (go.class in names(go)){
#            output_path_go <- paste0(go_dir, K, "_go_", go.class, "_", col_name, "_", cluster, ".tsv")
#            output_path_dotplot <- paste0(plots_dir, K, "_go_", go.class, "_", col_name, "_", cluster, ".pdf")
#            
#            go.df <- go[[go.class]]
#            go.df <- go.df %>%
#                  dplyr::mutate(
#                    Term = factor(Term, levels = rev(unique(Term))),
#                    signif = ifelse(FDR < FDR_thr, "Significant", "Not significant"),
#                    sizeVar = -log10(FDR)
#                  )
#
#            #p <- plotGO(go.df, FDR_thr)
#
#            write.table(x = go.df, file = output_path_go, col.names = T, row.names = F, sep = "\t", quote = F)
#            #ggplot2::ggsave(output_path_dotplot, p, width = 7, height = 7)
#          }
#        }
#
#    	  plot_autocorrelation <- plot_moran_autocorrelation(plot_data.s, K = K, plot_genes = F, plot_clusters = T)
#        output_path_plot_autocorrelation <- paste0(plots_dir, K, "_autocorrelation_annotated_", col_name, ".pdf")
#
#    	  ggplot2::ggsave(output_path_plot_autocorrelation, plot_autocorrelation, width = 7, height = 7)
#    }
#  )

local_moran_results <- lapply(X = (1:ncol(dataset))[-1],
                               FUN = function(idx){
                                     
                                     feature_name <- colnames(dataset)[idx]
                                                                           
                                     df <- data.frame(
                                      binID = dataset$binID,
                                      value = log1p(dataset[,idx])
                                     )

                                     localmoran_result <- localmoran(df$value, listw = listw)

                                     df$I <- localmoran_result[, "Ii"]
                                     df$pval <- localmoran_result[, "Pr(z != E(Ii))"]
                                     df$pval_fdr <- p.adjust(df$pval, "fdr")
                                     df$quadr <- attr(localmoran_result, "quadr")$mean                                     
                                     df$hotspot <- hotspot(obj = localmoran_result, Prname="Pr(z != E(Ii))", cutoff = 0.005, droplevels=FALSE)
                                     df$feature <- feature_name

                                     p <- visualize_local_moran(df, feature_name)
                                     
                                     ggplot2::ggsave(paste0(plots_dir, "/", "localmoran_", feature_name, ".png"), p, width = 14, height = 8)
                                     
                                     to.return <- list(df = df, p = p)

                                     return(to.return)

                               }
                              )
#}

bool_autocorrelation <- any(global_moran_results$FDR < FDR_thr)

writeLines(text = as.character(bool_autocorrelation), con = output_path_bool_autocorrelation, sep = "\n")
write.table(x = global_moran_results, file = output_path_univariate_table, sep = "\t", col.names = T, row.names = F, quote = F)
ggplot2::ggsave(output_path_plot_univariate, p, width = 7, height = 7)