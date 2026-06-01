source("workflow/scripts/helpers/find_neighbors.R")
source("workflow/scripts/helpers/smooth_something.R")

library(parallel)
library(dplyr)
library(tidyr)
library(ggplot2)

set.seed(123)
options(scipen=999)                                                               # unable scientific notation  

# args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)              # read args
# input_path_annotation               <- args$input_path_annotation
# output_path_annotation_smoothed	  <- args$output_path_annotation_smoothed
# B								      <- as.integer(args$B)						  # number of random shufflings
# C								      <- as.integer(args$C)						  # number of cores
# Ks_list				              <- seq(from = 1, to = 10)				      # K-neighbours to find for each side (e.g. K = 5 means 10-nearest-neighbours)


input_path_annotation                 <- ""
B								      <- 50
C								      <- 4
Ks_list				                  <- seq(from = 1, to = 10)				      # K-neighbours to find for each side (e.g. K = 5 means 10-nearest-neighbours)


`%>%` <- magrittr::`%>%`

annotation <- read.table(input_path_annotation, header = F, sep = "\t")
colnames(annotation) <- c("chr", "start", "end", "binID", "class")

smoothed_annotations <- parallel::mclapply(
						mc.cores = C, 
						X = 1:length(Ks_list),
						FUN = function(idx_out){

	K <- Ks_list[idx_out]
	annot_to_smooth <- annotation

	KNN <- parallel::mclapply(
		mc.cores = C, 
		X = 1:nrow(annot_to_smooth),
		FUN = function(idx_in){
				neighbors <- find_neighbors(annot_to_smooth, idx_in, K = K, return_self = T)
		}
	)

	annot_to_smooth$smoothed_class_final <- smooth_something(KNN, annot_to_smooth, "class", type = "char", C, B)$annotation
	smoothed_annot <- annot_to_smooth[,c("chr","start","end","binID","class","smoothed_class_final")]

	}
)

names(smoothed_annotations) <- as.character(Ks_list)

abundance_list <- list()

orig_counts <- annotation %>%
	dplyr::count(class) %>%
	dplyr::rename(abundance = n) %>%
	dplyr::mutate(K = "unsmoothed")

abundance_list[[1]] <- orig_counts

for(i in seq_along(smoothed_annotations)){

	Kname <- names(smoothed_annotations)[i]

	tmp <- smoothed_annotations[[i]] %>%
		dplyr::count(smoothed_class_final) %>%
		dplyr::rename(class = smoothed_class_final, abundance = n) %>%
		dplyr::mutate(K = paste0("K=", Kname))

	abundance_list[[length(abundance_list) + 1]] <- tmp
}

abundance_df <- dplyr::bind_rows(abundance_list)
abundance_df <- abundance_df %>% tidyr::complete(class, K, fill = list(abundance = 0))

abundance_df$K <- factor(abundance_df$K, levels = c("unsmoothed", paste0("K=", Ks_list)))

p <- ggplot(abundance_df, aes(x = K, y = abundance, fill = K)) +
	 geom_bar(stat = "identity") +
	 facet_wrap(~ class, scales = "free_y") +
	 theme_bw() +
	 theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(face = "bold"), legend.position = "none") +
	 labs(title = "Class abundance across smoothing levels", x = "Smoothing level", y = "Abundance")

print(p)



# write.table(
# 	x = annotation,
# 	file = output_path_annotation_smoothed,
# 	row.names = F,
# 	col.names = F,
# 	quote = F,
# 	sep = "\t"
# )