source("workflow/scripts/helpers/define_windows.R")                               # get helper functions
source("workflow/scripts/helpers/generate_backbone.R")

options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args


input_path_chrom_sizes         <- args$input_path_chrom_sizes
input_path_mask                <- args$input_path_mask
output_path_mask               <- args$output_path_mask
output_path_full               <- args$output_path_full
window_mode                    <- args$window_mode
window_size                    <- as.integer(args$window_size)
step                           <- as.integer(args$step)

chrom_sizes <- read.table(file = input_path_chrom_sizes,                          # read chrom sizes file
                          header = F,
                          sep = "\t",
                          col.names = c("chr", "size"))

chrom_sizes <- setNames(as.integer(chrom_sizes$size),                             # from df to named vector
                        as.character(chrom_sizes$chr))

windows <- define_windows(mode = window_mode,                                     # define Granges of chromosomes
                          step = step,
                          window_size = window_size,
                          chrom_sizes_canon = chrom_sizes)

backbone_full <- generate_backbone(windows)                                       # generate backbone

backbone_full$start <- backbone_full$start - 1                                    # set BED file coordinate (0-based, half-open)

if(file.exists(input_path_mask)){                                                 # if the mask file is provided, we use that one, otherwise, we use the full backbone as mask
  backbone_mask <- read.table(file = input_path_mask, header = F, sep = "\t")
} else {
  backbone_mask <- backbone_full
  stopifnot(identical(backbone_mask, backbone_full))
}

write.table(x = backbone_full,                                                    # write output
            file = output_path_full, 
            quote = F, sep = "\t", 
            row.names = F, col.names = F)

write.table(x = backbone_mask,                                                    # write output
            file = output_path_mask, 
            quote = F, sep = "\t", 
            row.names = F, col.names = F)