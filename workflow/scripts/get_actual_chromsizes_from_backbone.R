options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path                     <- args$input_path
output_path                    <- args$output_path

`%>%` <- magrittr::`%>%`

backbone <- read.table(input_path, header = F,  sep = "\t")
colnames(backbone) <- c("chr", "start", "end", "binID")

actual_chrom_size <- backbone %>% dplyr::group_by(chr) %>% dplyr::summarise(chrom_size = max(end))

write.table(file = output_path, x = actual_chrom_size, col.names = F, row.names = F, quote = F, sep = "\t")