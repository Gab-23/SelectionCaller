normATAC <- function(counts,minCount_acc,minSamples_acc,minCpm_acc,group,libSize,normalize){
	
	num_bins_before_filtering <- nrow(counts)
	
	keep_first_pass <- which(
				apply(X = counts[,-1], 
				      MARGIN = 1, 
				      minCount_acc = minCount_acc, 
				      minSamples_acc = minSamples_acc, 
				      function(x, minCount_acc, minSamples_acc){ 
				         sum(x > minCount_acc)}) > minSamples_acc)

	counts.f <- counts[keep_first_pass, ]
	num_bins_first_pass <- nrow(counts.f)

	if(!normalize){
	  return(counts.f)
	}

	counts.f <- DGEList(counts.f,group=group,lib.size=libSize)

	keep_second_pass <- rowSums(cpm(counts.f)>minCpm_acc) >= minSamples_acc

	counts.ff <- counts.f[keep_second_pass, , keep.lib.sizes=TRUE]

	num_bins_second_pass <- nrow(counts.ff)

	binIDs <- counts.ff$genes
	counts.ff.norm <- as.data.frame(cpm(counts.ff,normalized.lib.sizes=T)*mean(counts.ff$samples$lib.size)/100000)
	
	out.list <- list(binIDs = binIDs$binID, counts = counts.ff.norm)
	
	return(out.list)
}

# counts is my counts df
# minCount is an int
# minSamples is an int
# minCpm is an int
# group is a character vector
# libSize is a numerical vector
