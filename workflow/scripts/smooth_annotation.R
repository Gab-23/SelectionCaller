source("workflow/scripts/helpers/find_neighbors.R")				  # slight variation
source("workflow/scripts/helpers/smooth_something.R")

library(parallel)

set.seed(123)
options(scipen=999)                                                               # unable scientific notation  

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_annotation                 <- args$input_path_annotation
output_path_annotation_smoothed	      <- args$output_path_annotation_smoothed
K				                      <- as.integer(args$K) 					  # K-neighbours to find for each side (e.g. K = 5 means 10-nearest-neighbours)
B								      <- as.integer(args$B)						  # number of random shufflings
C								      <- as.integer(args$C)

`%>%` <- magrittr::`%>%`

annotation <- read.table(input_path_annotation, header = F, sep = "\t")
colnames(annotation) <- c("chr", "start", "end", "binID", "class")

KNN <- parallel::mclapply(
	   mc.cores = C, 
	   X = 1:nrow(annotation),
	   FUN = function(idx){
	           neighbors <- find_neighbors(annotation, idx, K = K, return_self = T)
	   }
	  )

# res <- do.call(
# 	cbind,
# 	parallel::mclapply(
# 	mc.cores = C,
# 	mc.set.seed = TRUE,
# 	X = seq(B), 
# 	FUN = function(b){
# 
# 		out <- rep(NA, nrow(annotation))
# 		
# 		shuffled_KNN <- KNN[sample.int(length(KNN), length(KNN))]
# 
# 		for (idx in seq(shuffled_KNN)) {
# 			shuffled_KNN.idxs <- shuffled_KNN[[idx]]
# 			unsmoothed_annot <- annotation$class[shuffled_KNN.idxs]
# 			most_common_class <- names(which.max(table(unsmoothed_annot)))
# 			out[shuffled_KNN.idxs] <- most_common_class
#     	}
# 
# 		return(out)
# 
# 		}
# 	)
# )
# 
# res <- as.data.frame(res)
# colnames(res) <- paste0("smoothed_class_", seq(B))
# annotation <- cbind(annotation, res)
# 
# stability <- apply(res, 1, function(x) {
#   max(table(x)) / length(x)
# })
# 
# annotation$smoothed_class_final <- lapply(X = seq(nrow(annotation)), 
# 									 	  FUN = function(idx){
# 									 		row <- unlist(annotation[idx, grep(pattern = "smoothed", x = colnames(annotation))])
# 									 		most_common_class <- names(which.max(table(row)))
# 									 		return(most_common_class)
# 									 	  }
# 										)

annotation$smoothed_class_final <- smooth_something(KNN, annotation, "class", type = "char", C, B)$annotation
annotation <- annotation[,c("chr","start","end","binID","class","smoothed_class_final")]

write.table(
	x = annotation,
	file = output_path_annotation_smoothed,
	row.names = F,
	col.names = F,
	quote = F,
	sep = "\t"
)