source("workflow/scripts/helpers/import_table.R")
source("workflow/scripts/helpers/process_dist_to_essential_genes_helper.R")

options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path              <- args$input_path
output_path             <- args$output_path

essential <- import_table(input_path, header = F)
colnames(essential) <- c("chr", "start", "end", "binID", "dist")

processed_essential <- process_dist_to_essential_genes(essential)
colnames(processed_essential) <- c("chr", "start", "end", "binID", "dist.to.closest.essential")

write.table(x = unique(processed_essential[,c("binID", "dist.to.closest.essential")]),    # write output
            file = output_path, 
            quote = F, 
            sep = "\t", 
            row.names = F, 
            col.names = T)