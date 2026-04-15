options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path                       <- args$input_path
window_size                      <- args$window_size
output_path                      <- args$output_path

backbone_annotated <- read.table(input_path, header = F, sep = "\t")
colnames(backbone_annotated) <- c("chr", "start", "end", "binID", "genes")

genes_on_bins <- lapply(X = 1:nrow(backbone_annotated), 
              FUN = function(idx){
              row <- backbone_annotated[idx,]; 
              genes <- strsplit(x = row$genes, split = ",")[[1]]; 
              return(genes)
           }
        )
names(genes_on_bins) <- backbone_annotated$binID

to.save <- list(genes_on_bins)
names(to.save) <- window_size

save(to.save, file = output_path)