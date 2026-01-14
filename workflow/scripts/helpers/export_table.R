export_table <- function(table, table_path, header = T, sep = "\t"){

    write.table(x = table,                                                         # write output
            file = table_path, 
            quote = F, 
            sep = sep, 
            row.names = F, 
            col.names = header)
}