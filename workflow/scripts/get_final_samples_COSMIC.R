options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args


input_path_samples                   <- args$input_path_samples
input_path_phenotype_ids             <- args$input_path_phenotype_ids
output_path_samples                  <- args$output_path_samples

phenotype_ids <- readLines(input_path_phenotype_ids)
samples <- read.delim(input_path_samples, header = T, sep = "\t")

samples_filtered_phenotype <- samples[samples$COSMIC_PHENOTYPE_ID %in% phenotype_ids,]

writeLines(text = samples_filtered_phenotype$COSMIC_SAMPLE_ID, sep = "\n", con = output_path_samples)