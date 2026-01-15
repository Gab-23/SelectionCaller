source("workflow/scripts/helpers/import_table.R")

options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_num_overlaps_protein_coding                         <- args$input_num_overlaps_protein_coding
input_num_overlaps_lncRNA                                 <- args$input_num_overlaps_lncRNA
input_num_overlaps_miRNA                                  <- args$input_num_overlaps_miRNA
input_num_overlaps_rRNA                                   <- args$input_num_overlaps_rRNA
output_num_overlaps_protein_coding                        <- args$output_num_overlaps_protein_coding
output_num_overlaps_lncRNA                                <- args$output_num_overlaps_lncRNA
output_num_overlaps_miRNA                                 <- args$output_num_overlaps_miRNA
output_num_overlaps_rRNA                                  <- args$output_num_overlaps_rRNA

gene_types <- c("protein_coding", "lncRNA", "miRNA", "rRNA")
input_overlaps <- list(input_num_overlaps_protein_coding, input_num_overlaps_lncRNA, input_num_overlaps_miRNA, input_num_overlaps_rRNA)
output_overlaps <- list(output_num_overlaps_protein_coding, output_num_overlaps_lncRNA, output_num_overlaps_miRNA, output_num_overlaps_rRNA)

for (idx in seq(length(input_overlaps))){

    gene_type <- gene_types[idx]
    input <- input_overlaps[[idx]]
    output <- output_overlaps[[idx]]
    
    table <- import_table(input, header = F)
    colnames(table) <- c("chr", "start", "end", "binID", paste0("num.overlapping.genes_", gene_type))
    
    write.table(x = table[,c("binID", paste0("num.overlapping.genes_", gene_type))],                                          # write output
            file = output, 
            quote = F, 
            sep = "\t", 
            row.names = F, 
            col.names = T)
}