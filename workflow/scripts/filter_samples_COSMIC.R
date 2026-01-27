options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args


input_path              <- args$input_path
output_path             <- args$output_path
sample_name_pattern     <- args$sample_name_pattern

samples <- read.delim(input_path, sep = "\t", header = T)

samples_filtered <- samples[which(startsWith(samples$SAMPLE_NAME, sample_name_pattern)),]

write.table(x = samples_filtered,                                                 # write output
            file = output_path, 
            quote = F, 
            sep = "\t", 
            row.names = F, 
            col.names = T)