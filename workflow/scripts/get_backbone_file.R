source("workflow/scripts/helpers/define_windows.R")                               # get helper functions
source("workflow/scripts/helpers/generate_backbone.R")

options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args


input_path       <- args$input_path
output_path      <- args$output_path
window_mode      <- args$window_mode
window_size      <- as.integer(args$window_size)
step             <- as.integer(args$step)

chrom_sizes <- read.table(file = input_path,                                      # read chrom sizes file
                          header = F,
                          sep = "\t",
                          col.names = c("chr", "size"))

chrom_sizes <- setNames(as.integer(chrom_sizes$size),                             # from df to named vector
                        as.character(chrom_sizes$chr))

windows <- define_windows(mode = window_mode,                                     # define Granges of chromosomes
                          step = step,
                          window_size = window_size,
                          chrom_sizes_canon = chrom_sizes)

backbone <- generate_backbone(windows)                                            # generate backbone

backbone$start <- backbone$start - 1                                              # set BED file coordinate (0-based, half-open)

write.table(x = backbone,                                                         # write output
            file = output_path, 
            quote = F, 
            sep = "\t", 
            row.names = F, 
            col.names = F)