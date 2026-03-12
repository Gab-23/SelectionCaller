source("workflow/scripts/helpers/donorm.R")
source("workflow/scripts/helpers/filterBins.R")

process_acc_vals <- function(counts, minCount_acc, minSamples_acc){
    AccValues <- counts
    rownames(AccValues) <- AccValues$binID
    AccValues$binID <- NULL
    AccValues.f <- AccValues #[, which(colSums(AccValues) >= 5000000)]
    AccValues.ff <- filterBins(AccValues.f, minCount=minCount_acc, minSamples=minSamples_acc)
    AccValues.tmm <- donorm(AccValues.ff)
    return(AccValues.tmm)
}