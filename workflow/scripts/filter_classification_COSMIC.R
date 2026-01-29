options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args


input_path                       <- args$input_path
output_path_phenotype_id         <- args$output_path_phenotype_id
COSMIC_parameters                <- args$COSMIC_parameters

COSMIC_parameters_parts <- strsplit(COSMIC_parameters, "_")[[1]]
site <- COSMIC_parameters_parts[1]
histology <- COSMIC_parameters_parts[2]

classification <- read.delim(input_path, sep = "\t", header = T)

classification_filtered <- classification[(classification$PRIMARY_SITE == site) & (classification$PRIMARY_HISTOLOGY == histology),]

#write.table(x = classification_filtered,                                          # write output
#            file = output_path_filtered, 
#            quote = F, 
#            sep = "\t", 
#            row.names = F, 
#            col.names = T)

writeLines(text = unique(classification_filtered$COSMIC_PHENOTYPE_ID), sep = "\n", con = output_path_phenotype_id)