source("workflow/scripts/helpers/normATAC.R")
source("workflow/scripts/helpers/pca_plot.R")


options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_dir                             <- args$input_dir
input_path_metadata                   <- args$input_path_metadata
minSamples_atacseq                    <- args$minSamples_atacseq                  # TODO: this can be automated in the snakefile, based on the metadata file
minCounts_atacseq                     <- args$minCounts_atacseq                   # TODO: this can be automated in the snakefile, based on the metadata file
minCpm_atacseq                        <- args$minCpm_atacseq                      # TODO: this can be automated in the snakefile, based on the metadata file
output_path_library_sizes             <- args$output_path_library_sizes
output_path_counts                    <- args$output_path_counts
output_path_PCA_condition             <- args$output_path_PCA_condition
output_path_PCA_library_size          <- args$output_path_PCA_library_size
output_path_PCA_sex                   <- args$output_path_PCA_sex
output_path_PCA_donor                 <- args$output_path_PCA_donor
output_path_PCA_study                 <- args$output_path_PCA_study
output_path_PCA_subtype               <- args$output_path_PCA_subtype 

files <- list.files(input_dir, recursive = T, full.names = T)
metadata <- read.table(input_path_metadata, header = T, sep = "\t")

sample_names <- unlist(lapply(X = files, FUN = function(x){strsplit(x, split = "/")[[1]][[5]]}))

cov_list <- lapply(X = seq(length(files)), 
                   FUN = function(idx){
                    file <- files[idx]
                    sample_name <- sample_names[idx]
                    tab <- read.table(file, header = F, sep = "\t")
                    colnames(tab) <- c("chr", "start", "end", "binID", sample_name)
                    to.return <- tab[,c("binID", sample_name)]
                    return(to.return)
                    })

counts <- Reduce(f = function(x,y){merge(x,y)}, 
                 x = cov_list)

library_sizes <- colSums(counts[,-1]) # in this case it makes sense to compute lib size like this because it's at the bin level. For peaks we need to use the .BAM files

group <- metadata$condition
names(group) <- metadata$SRA_id
group <- group[names(library_sizes)]

library_sizes_df <- data.frame(SRA_id = names(library_sizes),
                               library_size = library_sizes,
                               row.names = NULL)

metadata <- merge(metadata, library_sizes_df, by="SRA_id")

count.f.norm <- normATAC(counts = counts,
                         minSamples_acc = minSamples_atacseq,
                         minCount_acc = minCounts_atacseq,
                         minCpm_acc = minCpm_atacseq,
                         group = group,
                         libSize = library_sizes                       
                         )

pca <- prcomp(t(count.f.norm), scale. = TRUE)

pca_df <- data.frame(
  SRA_id = rownames(pca$x),
  PC1 = pca$x[,1],
  PC2 = pca$x[,2]
)

pca_df <- merge(pca_df, metadata, by = "SRA_id")

pca_plot_condition    <- pca_plot(pca_df, hue = "condition")          
pca_plot_library_size <- pca_plot(pca_df, hue = "library_size")
pca_plot_sex          <- pca_plot(pca_df, hue = "sex")                            # teniamo sex perchè potremmo volerlo per tts non BRCA
pca_plot_subtype      <- pca_plot(pca_df, hue = "subtype")      
pca_plot_donor        <- pca_plot(pca_df, hue = "donor_id", shape = "condition")
pca_plot_study        <- pca_plot(pca_df, hue = "study_id", shape = "condition")

write.table(x = count.f.norm,                                                     # write output
            file = output_path_counts, 
            quote = F, 
            sep = "\t", 
            row.names = F, 
            col.names = T)

write.table(x = library_sizes_df,
            file = output_path_library_sizes,
            quote = F,
            sep = "\t",
            row.names = F,
            col.names = T)

ggplot2::ggsave(output_path_PCA_condition, plot = pca_plot_condition, width = 6, height = 5)
ggplot2::ggsave(output_path_PCA_library_size, plot = pca_plot_library_size, width = 6, height = 5)
ggplot2::ggsave(output_path_PCA_sex, plot = pca_plot_sex, width = 6, height = 5)
ggplot2::ggsave(output_path_PCA_donor, plot = pca_plot_donor, width = 6, height = 5)
ggplot2::ggsave(output_path_PCA_study, plot = pca_plot_study, width = 6, height = 5)
ggplot2::ggsave(output_path_PCA_subtype, plot = pca_plot_subtype, width = 6, height = 5)