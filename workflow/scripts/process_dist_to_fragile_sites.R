source("workflow/scripts/helpers/import_table.R")
source("workflow/scripts/helpers/process_dist_to_fragile_sites_helper.R")

options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_dist_fragile_sites                                <- args$input_path_dist_fragile_sites
output_path_dist_fragile_sites                               <- args$output_path_dist_fragile_sites

table <- import_table(input_path_dist_fragile_sites, header = F)
colnames(table) <- c("chr","start","end", "binID", "dist")
    
processed_table <- process_dist_to_fragile_sites(table)
colnames(processed_table) <- c("chr", "start", "end", "binID", "dist.to.closest.FGS")

write.table(x = processed_table[,c("binID", "dist.to.closest.FGS")],                                          # write output
            file = output_path_dist_fragile_sites, 
            quote = F, 
            sep = "\t", 
            row.names = F, 
            col.names = T)
