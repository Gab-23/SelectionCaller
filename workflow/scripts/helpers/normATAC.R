normATAC <- function(counts,minCount_acc,minSamples_acc,minCpm_acc,group,libSize){
	num_bins_before_filtering <- nrow(counts)
	
	o <- counts[which(
				  apply(X = counts[,-1], 
					    MARGIN = 1, 
					    minCount_acc = minCount_acc, 
					    minSamples_acc = minSamples_acc, 
					    function(x, minCount_acc, minSamples_acc){ 
					  	  sum(x > minCount_acc)
					  	}) > minSamples_acc), ]

	num_bins_second_step <- nrow(o)
	
	o <- DGEList(o,group=group,lib.size=libSize)
	keep <- rowSums(cpm(o)>minCpm_acc) >= minSamples_acc
	o <- o[keep, , keep.lib.sizes=TRUE]

	num_bins_third_step <- nrow(o)

	o <- as.data.frame(cpm(o,normalized.lib.sizes=T)*mean(o$samples$lib.size)/10000)
}

# counts is my counts df
# minCount is an int
# minSamples is an int
# minCpm is an int
# group is a character vector
# libSize is a numerical vector