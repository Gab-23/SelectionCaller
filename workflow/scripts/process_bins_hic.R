options(scipen=999)                                                               # unable scientific notation

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_backbone                     <- args$input_path_backbone
input_path_bins                         <- args$input_path_bins
output_path                             <- args$output_path

`%>%` <- magrittr::`%>%`


backbone <- read.table(file = input_path_backbone, header = F, sep = "\t")
bins <- read.table(file = input_path_bins, header = F, sep = "\t")

colnames(backbone) <- c("chr", "start", "end", "binID")
colnames(bins) <- c("chr", "start", "end", "unk")

bins_with_binID <- merge(x = backbone, y = bins, by = c("chr", "start", "end"))
bins_with_binID$unk <- NULL

write.table(x = bins_with_binID, file = output_path, quote = F, col.names = F, row.names = F, sep = "\t")
