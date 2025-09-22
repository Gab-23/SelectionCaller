source("workflow/scripts/define_windows.R")                                       # get helper functions
source("workflow/scripts/generate_backbone.R")

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args


input_file       <- args$input_file
output_file      <- args$output_file
window_mode      <- args$window_mode
window_size      <- as.integer(args$window_size)
step             <- as.integer(args$step)

chrom_canon <- paste0("chr",1:22)                                                 # define desired chromosomes

chrom_sizes <- read.table(file = input_file,
                          header = F,
                          sep = "\t",
                          col.names = c("chr", "size"))

chrom_sizes_canon <- chrom_sizes[chrom_sizes$chr %in% chrom_canon,]
chrom_sizes_canon <- setNames(as.integer(chrom_sizes_canon$size), 
                              as.character(chrom_sizes_canon$chr))

windows <- define_windows(mode = window_mode,
                          step = step,
                          window_size = window_size,
                          chrom_sizes_canon = chrom_sizes_canon)

backbone <- generate_backbone(windows)

backbone$start <- backbone$start - 1

write.table(x = backbone, 
            file = output_file, 
            quote = F, 
            sep = "\t", 
            row.names = F, 
            col.names = F)