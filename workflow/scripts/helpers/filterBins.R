filterGenes <- function(counts, minCount=20, minSamples=10){
	nb1 <- nrow(counts)
	counts <- counts[which(apply(counts,1,minCount=minCount,minSamples=minSamples,function(x,minCount,minSamples){ sum(x>minCount) })>minSamples),, drop = FALSE]
	message(paste(nb1-nrow(counts),"rows out of",nb1,"were removed. [COUNTS]"))
	return(counts)
}