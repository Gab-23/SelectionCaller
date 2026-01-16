source("workflow/scripts/helpers/import_table.R")

options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path                <- args$input_path
output_path               <- args$output_path

gc_content.df <- import_table(input_path, header = F)
colnames(gc_content.df) <- c("chr", "start", "end", "binID", "GC_content")

write.table(x = gc_content.df[,c("binID", "GC_content")],                                                     # write output
            file = output_path, 
            quote = F, 
            sep = "\t", 
            row.names = F, 
            col.names = T)