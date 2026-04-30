# this assumes that the dataset is chromosome sorted
source("workflow/scripts/helpers/do_simple_boxplots.R")
source("workflow/scripts/helpers/find_neighbors.R")				  # slight variation
source("workflow/scripts/helpers/smooth_something.R")
source("workflow/scripts/helpers/boxplot_0s_non0s.R")

library(parallel)

set.seed(123)
options(scipen=999)                                                               # unable scientific notation

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

plots_dir                                            <- args$plots_dir
input_path_dataset                                   <- args$input_path_dataset
output_path_boxplots_0s_non0s_muts                   <- args$output_path_boxplots_0s_non0s_muts
output_path_boxplots_0s_non0s_chromstates            <- args$output_path_boxplots_0s_non0s_chromstates
output_path_boxplots_0s_non0s_muts_smooth            <- args$output_path_boxplots_0s_non0s_muts_smooth
output_path_boxplots_0s_non0s_chromstates_smooth     <- args$output_path_boxplots_0s_non0s_chromstates_smooth

K				                                     <- as.integer(args$K) 					  # K-neighbours to find for each side (e.g. K = 5 means 10-nearest-neighbours)
B								                     <- as.integer(args$B)						  # number of random shufflings
C								                     <- as.integer(args$C)
pattern1                                             <- "^length_(non)?coding"
pattern2                                             <- "^[0-9]+_[aA-zZ]+"

`%>%` <- magrittr::`%>%`

dataset <- read.delim(input_path_dataset, check.names = F, header = T, sep = "\t")
colnames(dataset) <- gsub(pattern = "/", replacement = "_", x = colnames(dataset)) # clean colnames of chromatin states. 20_ZNF/Rpts -> 20_ZNF_Rpts

cols_to_smooth <- grep(pattern = pattern1, x = colnames(dataset), value = T)
cols_to_smooth <- c(cols_to_smooth, grep(pattern = pattern2, x = colnames(dataset), value = T))

unique_chromosomes <- paste0(1:22) 

for(chr in unique_chromosomes){

    dataset.f <- dataset[startsWith(x = dataset$binID, prefix = paste0(chr, "_")),]

    list_of_distributions <- list()
    # list_of_sd <- list()

    for (col in cols_to_smooth){
            
        col_of_interest <- dataset.f[,col]
        col_of_interest <- ifelse(col_of_interest == 0, 0, 1)
        col_of_interest <- ifelse(is.na(col_of_interest), 2, col_of_interest)

        zeros_idxs <- which(col_of_interest == 0)

        stretch_of_similar_values <- cumsum(c(1, as.numeric(diff(col_of_interest) != 0)))
        stretch_of_similar_values.f <- stretch_of_similar_values[zeros_idxs]

        tab <- table(stretch_of_similar_values.f)
        names(tab) <- NULL

        distribution_of_zeros <- as.numeric(sort(tab))

        list_of_distributions[[col]] <- distribution_of_zeros
        # list_of_sd[[col]] <- ifelse(is.na(sd(distribution_of_zeros)), 0, sd(distribution_of_zeros))

        }

        output_path_boxplot <- paste0(plots_dir, "/", chr, "_contiguous_zero_distribution.png")

        p <- do_simple_boxplots(list_of_distributions,
                                x_axis = "stretches of 0s", 
                                title = paste0("chromosome", chr)
                            )

        ggplot2::ggsave(output_path_boxplot, p, width = 7, height = 7)

    }

KNN <- parallel::mclapply(
	   mc.cores = C, 
	   X = 1:nrow(dataset),
	   FUN = function(idx){
	           neighbors <- find_neighbors(dataset, idx, K = K, return_self = T)
	   }
	  )

smoothed_cols <- do.call(cbind,
        parallel::mclapply(
            mc.cores = C,
            X = cols_to_smooth,
            FUN = function(col){
                smoothed <- smooth_something(KNN, dataset, col, type = "num", C, B)$smoothed
                smoothed.df <- data.frame(val = smoothed)
                colnames(smoothed.df) <- paste0(col, "_smoothed")
                print(col)
                return(smoothed.df)
            }
        ))

dataset <- cbind(dataset, smoothed_cols)

#for (col in cols_to_smooth){
#    dataset[, paste0(col, "_smoothed")] <- smooth_something(KNN, dataset, col, type = "num", C, B)$smoothed
#}

boxplots_0s_non0s_muts <- boxplots_0s_non0s(dataset[,c("binID", "logFC", cols_to_smooth)], pattern1)
boxplots_0s_non0s_muts_smooth <- boxplots_0s_non0s(dataset, paste0(pattern1, ".*smoothed$"))

boxplots_0s_non0s_chromstates <- boxplots_0s_non0s(dataset[,c("binID", "logFC", cols_to_smooth)], pattern2)
boxplots_0s_non0s_chromstates_smooth <- boxplots_0s_non0s(dataset, paste0(pattern2, ".*smoothed$"))

ggplot2::ggsave(output_path_boxplots_0s_non0s_muts, boxplots_0s_non0s_muts, width = 20, height = 20)
ggplot2::ggsave(output_path_boxplots_0s_non0s_muts_smooth, boxplots_0s_non0s_muts_smooth, width = 20, height = 20)
ggplot2::ggsave(output_path_boxplots_0s_non0s_chromstates, boxplots_0s_non0s_chromstates, width = 20, height = 20)
ggplot2::ggsave(output_path_boxplots_0s_non0s_chromstates_smooth, boxplots_0s_non0s_chromstates_smooth, width = 20, height = 20)
