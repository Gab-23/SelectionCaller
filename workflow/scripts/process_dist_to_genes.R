source("workflow/scripts/helpers/import_table.R")
source("workflow/scripts/helpers/process_dist_to_genes_helper.R")

options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_dist_protein_coding                                 <- args$input_dist_protein_coding
input_dist_lncRNA                                         <- args$input_dist_lncRNA
input_dist_miRNA                                          <- args$input_dist_miRNA
input_dist_rRNA                                           <- args$input_dist_rRNA   
output_dist_protein_coding                                <- args$output_dist_protein_coding
output_dist_lncRNA                                        <- args$output_dist_lncRNA
output_dist_miRNA                                         <- args$output_dist_miRNA
output_dist_rRNA                                          <- args$output_dist_rRNA

gene_types <- c("protein_coding", "lncRNA", "miRNA", "rRNA")
input_dist <- list(input_dist_protein_coding, input_dist_lncRNA, input_dist_miRNA, input_dist_rRNA)
output_dist <- list(output_dist_protein_coding, output_dist_lncRNA, output_dist_miRNA, output_dist_rRNA)

for (idx in seq(length(input_dist))){

    gene_type <- gene_types[idx]
    input <- input_dist[[idx]]
    output <- output_dist[[idx]]
    
    table <- import_table(input, header = F)
    colnames(table) <- c("chr","start","end", "binID", "ID", "gene_name", "type", "dist")
    
    processed_table <- process_dist_to_genes(table)
    colnames(processed_table) <- c("chr","start","end", "binID", "ID", "gene_name", "type", paste0("dist.to.closest.gene_", gene_type))
    
    write.table(x = processed_table[,c("binID", paste0("dist.to.closest.gene_", gene_type))],                                          # write output
            file = output, 
            quote = F, 
            sep = "\t", 
            row.names = F, 
            col.names = T)
}