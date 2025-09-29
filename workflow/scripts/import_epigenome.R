import_epigenome <- function(epigenome_path){

    epigenome <- read.table(file = epigenome_path,
                            header = TRUE, 
                            sep = "\t", 
                            quote = "",
                            check.names = FALSE)

    return(epigenome)
}