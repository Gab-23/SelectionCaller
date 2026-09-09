options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path             <- args$input_path
output_path            <- args$output_path

chrom_sizes_canon <- read.table(input_path)

starts <- data.frame(chr = chrom_sizes_canon$V1, start = rep(0, nrow(chrom_sizes_canon)), end = rep(1, nrow(chrom_sizes_canon)))

ends <- chrom_sizes_canon
ends$V3 <- ends$V2
ends$V2 <- ends$V2 - 1
colnames(ends) <- c("chr", "start", "end")

telomeres <- rbind(starts, ends)

write.table(x = telomeres[order(telomeres$chr),],                                 # write output
            file = output_path, 
            quote = F, 
            sep = "\t", 
            row.names = F, 
            col.names = F)