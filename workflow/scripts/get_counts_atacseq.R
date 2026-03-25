source("workflow/scripts/helpers/normATAC.R")
source("workflow/scripts/helpers/svaBatchCor.R")
source("workflow/scripts/helpers/pca_plot.R")

library(edgeR)
library(sva)

options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_dir                             <- args$input_dir
input_path_metadata                   <- args$input_path_metadata
input_path_samples_to_exclude         <- args$input_path_samples_to_exclude

minSamples_atacseq                    <- args$minSamples_atacseq                  # TODO: this can be automated in the snakefile, based on the metadata file
minCounts_atacseq                     <- args$minCounts_atacseq                   # TODO: this can be automated in the snakefile, based on the metadata file
minCpm_atacseq                        <- args$minCpm_atacseq                      # TODO: this can be automated in the snakefile, based on the metadata file
SVs				      <- args$SVs

output_path_library_sizes             <- args$output_path_library_sizes

output_path_counts                    <- args$output_path_counts
output_path_counts_sva                <- args$output_path_counts_sva

output_path_PCA_condition             <- args$output_path_PCA_condition
output_path_PCA_library_size          <- args$output_path_PCA_library_size
output_path_PCA_sex                   <- args$output_path_PCA_sex
output_path_PCA_donor                 <- args$output_path_PCA_donor
output_path_PCA_study                 <- args$output_path_PCA_study
output_path_PCA_subtype               <- args$output_path_PCA_subtype 

output_path_PCA_sva_condition          <- args$output_path_PCA_sva_condition
output_path_PCA_sva_library_size       <- args$output_path_PCA_sva_library_size
output_path_PCA_sva_sex                <- args$output_path_PCA_sva_sex
output_path_PCA_sva_donor              <- args$output_path_PCA_sva_donor
output_path_PCA_sva_study              <- args$output_path_PCA_sva_study
output_path_PCA_sva_subtype            <- args$output_path_PCA_sva_subtype 

files <- list.files(input_dir, recursive = T, full.names = T)
metadata <- read.table(input_path_metadata, header = T, sep = "\t")
samples_to_exclude <- readLines(input_path_samples_to_exclude)

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

counts <- counts[,!(colnames(counts) %in% samples_to_exclude)]

library_sizes <- colSums(counts[,-1]) # in this case it makes sense to compute lib size like this because it's at the bin level. For peaks we need to use the .BAM files

library_sizes_df <- data.frame(SRA_id = names(library_sizes),
                               library_size = library_sizes,
                               row.names = NULL)

group <- metadata$condition
names(group) <- metadata$SRA_id
group <- group[names(library_sizes)]

metadata <- merge(metadata, library_sizes_df, by="SRA_id")

normATAC.out <- normATAC(counts = counts,
                         minSamples_acc = minSamples_atacseq,
                         minCount_acc = minCounts_atacseq,
                         minCpm_acc = minCpm_atacseq,
                         group = group,
                         libSize = library_sizes                       
                         )

count.f.norm <- normATAC.out$counts
binIDs <- normATAC.out$binIDs
count.f.norm$binID <- binIDs
count.f.norm <- count.f.norm[, c("binID", setdiff(names(count.f.norm), "binID"))]

########## PRE SVA

pca <- prcomp(t(count.f.norm[,-1]), scale. = TRUE)

pca_df <- data.frame(
  SRA_id = rownames(pca$x),
  PC1 = pca$x[,1],
  PC2 = pca$x[,2]
)

pca_df <- merge(pca_df, metadata, by = "SRA_id")

pca_plot_condition    <- pca_plot(pca, pca_df, hue = "condition")          
pca_plot_library_size <- pca_plot(pca, pca_df, hue = "library_size")
pca_plot_sex          <- pca_plot(pca, pca_df, hue = "sex")                            # teniamo sex perchè potremmo volerlo per tts non BRCA
pca_plot_subtype      <- pca_plot(pca, pca_df, hue = "subtype")      
pca_plot_donor        <- pca_plot(pca, pca_df, hue = "donor_id", shape = "condition")
pca_plot_study        <- pca_plot(pca, pca_df, hue = "study_id", shape = "condition")

############ POST SVA

# definisci il batch che credi STA
study <- metadata$study_id

# definisci la variabile INDIE che vuoi testare davvero (in questo caso il genotipo)
condition <- metadata$condition

# model matrix senza covariate così vediamo quanta merda c'è
mm1 <- model.matrix(~condition + study)

# conta il numero di surrogate variables
n.sv <- num.sv(count.f.norm[,-1], mm1, method="leek")

# model matrix con solo le covariate
mm0 <- model.matrix(~study)

# salvo da parte i valori di SV
sv <- sva(as.matrix(count.f.norm[,-1]),mm1,mm0)$sv

## generi la model matrix magggica
mm <- model.matrix(~sv+condition)

# correggi le conte con la funzia di PL (nota sempre il 2)
count.f.norm.sva <- svaBatchCor(count.f.norm[,-1], mm1, mm0, SVs)$corrected

pca_sva <- prcomp(t(count.f.norm.sva), scale. = TRUE)

pca_sva_df <- data.frame(
  SRA_id = rownames(pca_sva$x),
  PC1 = pca_sva$x[,1],
  PC2 = pca_sva$x[,2]
)

pca_sva_df <- merge(pca_sva_df, metadata, by = "SRA_id")

pca_sva_plot_condition    <- pca_plot(pca_sva, pca_sva_df, hue = "condition")          
pca_sva_plot_library_size <- pca_plot(pca_sva, pca_sva_df, hue = "library_size")
pca_sva_plot_sex          <- pca_plot(pca_sva, pca_sva_df, hue = "sex")                            # teniamo sex perchè potremmo volerlo per tts non BRCA
pca_sva_plot_subtype      <- pca_plot(pca_sva, pca_sva_df, hue = "subtype")      
pca_sva_plot_donor        <- pca_plot(pca_sva, pca_sva_df, hue = "donor_id", shape = "condition")
pca_sva_plot_study        <- pca_plot(pca_sva, pca_sva_df, hue = "study_id", shape = "condition")

############

write.table(x = count.f.norm,                                                     # write output
            file = output_path_counts, 
            quote = F, 
            sep = "\t", 
            row.names = F, 
            col.names = T)

write.table(x = count.f.norm.sva,                                                 # write output
            file = output_path_counts_sva, 
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

ggplot2::ggsave(output_path_PCA_sva_condition, plot = pca_sva_plot_condition, width = 6, height = 5)
ggplot2::ggsave(output_path_PCA_sva_library_size, plot = pca_sva_plot_library_size, width = 6, height = 5)
ggplot2::ggsave(output_path_PCA_sva_sex, plot = pca_sva_plot_sex, width = 6, height = 5)
ggplot2::ggsave(output_path_PCA_sva_donor, plot = pca_sva_plot_donor, width = 6, height = 5)
ggplot2::ggsave(output_path_PCA_sva_study, plot = pca_sva_plot_study, width = 6, height = 5)
ggplot2::ggsave(output_path_PCA_sva_subtype, plot = pca_sva_plot_subtype, width = 6, height = 5)
