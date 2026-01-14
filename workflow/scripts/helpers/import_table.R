import_table <- function(table_path, header = TRUE, sep = "\t"){

    table <- read.table(file = table_path,
                            header = header, 
                            sep = sep, 
                            quote = "",
                            check.names = FALSE)

    return(table)
}