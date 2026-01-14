source("workflow/scripts/helpers/import_table.R")
source("workflow/scripts/helpers/process_dist_to_cancer_genes_helper.R")

options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_TSG              <- args$input_path_TSG
input_path_OG               <- args$input_path_OG
output_path_TSG             <- args$output_path_TSG
output_path_OG              <- args$output_path_OG

TSG <- import_table(input_path_TSG, header = F)
colnames(TSG) <- c("chr", "start", "end", "binID", "gene", "role", "dist")

OG <- import_table(input_path_OG, header = F) 
colnames(OG) <- c("chr", "start", "end", "binID", "gene", "role", "dist")

processed_TSG <- process_dist_to_cancer_genes(TSG)
colnames(processed_TSG) <- c("chr", "start", "end", "binID", "gene", "role", "dist.to.closest.TSG")

processed_OG <- process_dist_to_cancer_genes(OG)
colnames(processed_OG) <- c("chr", "start", "end", "binID", "gene", "role", "dist.to.closest.OG")

write.table(x = processed_TSG[,c("binID", "dist.to.closest.TSG")],                                                 # write output
            file = output_path_TSG, 
            quote = F, 
            sep = "\t", 
            row.names = F, 
            col.names = T)

write.table(x = processed_OG[,c("binID", "dist.to.closest.OG")],                                                 # write output
            file = output_path_OG, 
            quote = F, 
            sep = "\t", 
            row.names = F, 
            col.names = T)