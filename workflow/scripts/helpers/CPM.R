CPM <- function(df, library_size){
    #CPM = readsMappedToGene * 1/totalNumReads * 10^6
    df$CPM <- ((df$coverage) * (1 / library_size) * (10**6))
    return(df)
}