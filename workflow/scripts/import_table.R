import_table <- function(table_path){

    table <- read.table(file = table_path,
                            header = TRUE, 
                            sep = "\t", 
                            quote = "",
                            check.names = FALSE)

    return(table)
}