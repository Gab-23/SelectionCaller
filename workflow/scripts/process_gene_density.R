options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_essential             <- args$input_path_essential
input_path_OGs                   <- args$input_path_OGs
input_path_TSGs                  <- args$input_path_TSGs
output_path_essential            <- args$output_path_essential
output_path_OGs                  <- args$output_path_OGs
output_path_TSGs                 <- args$output_path_TSGs

essential <- read.table(input_path_essential, sep = "\t", header = F)
colnames(essential) <- c("chr", "start", "end", "binID", "genes", "target_gene")

OGs <- read.table(input_path_OGs, sep = "\t", header = F)
colnames(OGs) <- c("chr", "start", "end", "binID", "genes", "target_gene")

TSGs <- read.table(input_path_OGs, sep = "\t", header = F)
colnames(TSGs) <- c("chr", "start", "end", "binID", "genes", "target_gene")

data_list <- list(essential = essential, OGs = OGs, TSGs = TSGs)
output_paths_list <- list(essential = output_path_essential, OGs = output_path_OGs, TSGs = output_path_TSGs)

for (data_type in names(data_list)){
    
    data <- data_list[[data_type]]
    output_path <- output_paths_list[[data_type]]

    data$num_genes <- unlist(lapply(X = data$genes, FUN = function(x){length(strsplit(x, split=",")[[1]])}))
    data$num_target_genes <- unlist(lapply(X = data$target_gene, FUN = function(x){length(strsplit(x, split=",")[[1]])}))
    data[, paste0("density_", data_type)] <- data$num_target_genes / data$num_genes

    write.table(x = data, file = output_path, col.names = T, row.names = T, quote = F, sep = "\t")

}