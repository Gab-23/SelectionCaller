source("workflow/scripts/helpers/import_table.R")

options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_num_overlaps_fragile_sites                         <- args$input_path_num_overlaps_fragile_sites
output_path_num_overlaps_fragile_sites                        <- args$output_path_num_overlaps_fragile_sites

table <- import_table(input_path_num_overlaps_fragile_sites, header = F)
colnames(table) <- c("chr", "start", "end", "binID", "num.overlapping.FGS")
    
write.table(x = table[,c("binID", "num.overlapping.FGS")],                                          # write output
            file = output_path_num_overlaps_fragile_sites, 
            quote = F, 
            sep = "\t", 
            row.names = F, 
            col.names = T)
