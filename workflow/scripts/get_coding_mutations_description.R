options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_coords              <- args$input_path_coords
input_path_muts                <- args$input_path_muts
output_path                    <- args$output_path

coords <- read.delim("resources/COSMIC/coordinates/coding/coordinates.sorted.bed", header = F)
muts <- read.delim("resources/COSMIC/mutations/coding/mutations.filtered.bed", header = T)

colnames(coords) <- c("chr", "start", "end", "mutID")

muts_sub <- unique(muts[,c("CHROMOSOME","GENOME_START","GENOME_STOP","GENOMIC_MUTATION_ID","MUTATION_DESCRIPTION")])

merged <- merge(coords, muts_sub, by.x = "mutID", by.y = "GENOMIC_MUTATION_ID")

output.df <- merged[,c("mutID", "chr", "start", "end", "MUTATION_DESCRIPTION")]
colnames(output.df) <- c("mutID", "chr", "start", "end", "mut_types")

output.df <- output.df[,c("chr", "start", "end", "mutID", "mut_types")]

write.table(x = output.df,                                                        # write output
            file = output_path, 
            quote = F, 
            sep = "\t", 
            row.names = F, 
            col.names = F)