options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path              <- args$input_path
output_path             <- args$output_path
essentiality_cutoff     <- args$essentiality_cutoff
num_assays_cutoff       <- args$num_assays_cutoff

crispr <- read.delim(input_path, sep = ",")

col_names <- gsub(pattern = "\\.\\..*", replacement = "", x = colnames(crispr))
col_names <- gsub(pattern = "\\.", replacement = "-", x = col_names)
col_names[1] <- "ModelID"

colnames(crispr) <- col_names
crispr[is.na(crispr)] <- 0

essential_genes_per_assay <- unlist(apply(X = crispr[,-1], 
                                   MARGIN = 1, 
                                   FUN = function(row){
                                       names <- names(row)
                                       row <- as.numeric(row)
                                       names(row) <- names
                                       return(names(row[row>essentiality_cutoff]))
                                       }
                                   ))

num_essential_assays_per_gene <- table(essential_genes_per_assay)

essential_genes <- names(num_essential_assays_per_gene[num_essential_assays_per_gene > num_assays_cutoff])

writeLines(text = essential_genes, con = output_path, sep = "\n")
