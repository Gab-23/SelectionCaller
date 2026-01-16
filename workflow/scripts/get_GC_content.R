# taken from https://www.biostars.org/p/478444/

options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path                <- args$input_path
output_path               <- args$output_path
assembly_code             <- args$assembly_code

genome_name <- paste0("BSgenome.Hsapiens.UCSC.", assembly_code)
genome <- BSgenome::getBSgenome(genome_name)

backbone.gr <- rtracklayer::import.bed(input_path)

gc_content <- Repitools::gcContentCalc(backbone.gr , organism = genome, verbose = F)

GenomicRanges::mcols(backbone.gr)$GC_content <- gc_content

gc_content.df <- as.data.frame(backbone.gr)
gc_content.df$width <- NULL
gc_content.df$strand <- NULL

gc_content.df$start <- gc_content.df$start - 1                                     # set BED file coordinate (0-based, half-open)

write.table(x = gc_content.df,                                                     # write output
            file = output_path, 
            quote = F, 
            sep = "\t", 
            row.names = F, 
            col.names = F)