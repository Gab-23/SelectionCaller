# idea: generate metric for each dataset and evaluate

source("workflow/scripts/helpers/find_neighbors.R")
source("workflow/scripts/helpers/connect_KNN.R")
source("workflow/scripts/helpers/connect_KNN_igraph.R")

options(scipen=999)                                                               # unable scientific notation           

library(FNN)
library(ComplexHeatmap)

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_dataset             <- args$input_path_dataset
input_path_backbone            <- args$input_path_backbone
C                              <- as.integer(args$C)

dataset <- read.delim(input_path_dataset, check.names = F, header = T, sep = "\t")
colnames(dataset) <- gsub(pattern = "/", replacement = "_", x = colnames(dataset))

backbone <- read.table(input_path_backbone, header = F, sep = "\t")
colnames(backbone) <- c("chr", "start", "end", "binID")

unique_chromosomes <- 1:22
K_list <- 1:10

KNNs_list <- parallel::mclapply(
    mc.cores = C, 
    X = K_list,
    FUN = function(K){

    KNN <- parallel::mclapply(
        mc.cores = C,
        X = unique_chromosomes,
        FUN = function(chr){
            
        backbone.f <- backbone[startsWith(prefix = paste0(chr, "_"), x = backbone$binID),]
        dataset.f <- dataset[dataset$binID %in% backbone.f$binID,]
        backbone.f <- backbone.f[backbone.f$binID %in% dataset.f$binID,]

        all(backbone.f$binID == dataset.f$binID)

        KNN.chr <- parallel::mclapply(
                    mc.cores = C,
                    X = 1:nrow(backbone.f), 
                    FUN = function(idx){

                        row <- dataset.f[idx,]
                        
                        neighbors.spatial <- find_neighbors(backbone.f, idx, K, F)
                        
                        neighbors.data <- FNN::get.knnx(dataset.f, row, ((K*2)+1), "kd_tree")$nn.index[1,]
                        neighbors.data <- neighbors.data[neighbors.data != idx]
                        # if (data_mode == "corr"){}else{}
                        # K*2 is because find_neighbors looks for K upstream and K downstream
                        # here we do need to use FNN::get.knnx(dataset, row, ((K*2)+1), "kd_tree")$nn.index[1,] since we do not have the correlation space
                        
                        neighbors.intersect <- intersect(neighbors.spatial, neighbors.data)
                        if (length(neighbors.intersect) > 0){
                            return(neighbors.intersect)
                        } else {
                            return(NULL)
                        }
                    })

        names(KNN.chr) <- backbone.f$binID
        
        KNN.chr.connected <- connect_KNN_igraph(KNN.chr)
        names(KNN.chr.connected) <- paste0(chr, "_", names(KNN.chr.connected))

        return(KNN.chr.connected)
        
        })

    names(KNN) <- unique_chromosomes
    KNN <- do.call(c, KNN) # flatten

    freq.closest.similar.bins <- lapply(X = KNN, 
                                        FUN = function(x){
                                            to.return <- (length(x)/K)
                                            return(to.return)
                                        })
    
    avg.freq.closest.similar.bins <- mean(unlist(freq.closest.similar.bins))

    to.return <- list(KNN = KNN, statistic = statistic)

    return(to.return)

    })

