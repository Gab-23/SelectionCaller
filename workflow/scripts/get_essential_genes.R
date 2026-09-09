source("workflow/scripts/helpers/parse_depmap_colnames.R")

options(scipen=999)                                                               # unable scientific notation

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path              <- args$input_path
scores_path             <- args$scores_path
output_path             <- args$output_path
essentiality_cutoff     <- args$essentiality_cutoff
num_assays_cutoff       <- args$num_assays_cutoff

crispr <- read.delim(input_path, sep = ",", header = T)
colnames(crispr) <- parse_depmap_gene_symbols(colnames(crispr))
crispr[is.na(crispr)] <- 0

essential_genes_per_assay <- unlist(apply(X = crispr[,-1],
                                   MARGIN = 1,
                                   FUN = function(row){
                                       names <- names(row)
                                       row <- as.numeric(row)
                                       names(row) <- names
                                       return(names(row[row > essentiality_cutoff]))
                                       }
                                   ))

num_essential_assays_per_gene <- table(essential_genes_per_assay)
essential_genes <- names(num_essential_assays_per_gene[num_essential_assays_per_gene > num_assays_cutoff])

scores_header <- colnames(read.delim(scores_path, sep = ",", header = TRUE, nrows = 0))
scores_genes  <- parse_depmap_gene_symbols(scores_header)

dropped_genes <- setdiff(essential_genes, scores_genes)
if (length(dropped_genes) > 0) {
    message(sprintf(
        "get_essential_genes.R: dropping %d gene(s) absent from scores file (likely readthrough/fusion genes): %s",
        length(dropped_genes), paste(dropped_genes, collapse = ", ")
    ))
}
essential_genes <- intersect(essential_genes, scores_genes)

writeLines(text = essential_genes, con = output_path, sep = "\n")