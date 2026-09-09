parse_depmap_gene_symbols <- function(colnames_vec) {
    cn <- gsub(pattern = "\\.\\..*", replacement = "", x = colnames_vec)
    cn <- gsub(pattern = "\\.", replacement = "-", x = cn)
    cn[1] <- "ModelID"
    return(cn)
}