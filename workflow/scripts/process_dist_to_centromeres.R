source("workflow/scripts/helpers/import_table.R")

options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path                <- args$input_path
output_path               <- args$output_path

dist_to_centromeres <- import_table(input_path, header = F)
colnames(dist_to_centromeres) <- c("chr", "start", "end", "binID", "dist_to_centromere")

write.table(x = dist_to_centromeres[,c("binID", "dist_to_centromere")],                   # write output
            file = output_path, 
            quote = F, 
            sep = "\t", 
            row.names = F, 
            col.names = T)