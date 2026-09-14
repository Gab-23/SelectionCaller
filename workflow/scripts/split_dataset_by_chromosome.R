options(scipen=999)                                                               # unable scientific notation  

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path                          <- args$input_path
chrom                               <- args$chrom 
output_path                         <- args$output_path

pattern <- paste0(chrom, "_")

dataset <- read.table(input_path, header = T, sep = "\t")
dataset.f <- dataset[startsWith(dataset$binID, pattern), , drop = FALSE]

write.table(x = dataset.f, file = output_path, sep = "\t", row.names = F, col.names = T, quote = F)