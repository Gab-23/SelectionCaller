options(scipen=999)                                                              # unable scientific notation     

library(parallel)
library(stringr)

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)               # read args

# backbone_table                       <- args$backbone_table                      # "features_df.tsv"
# window_size                          <- args$window_size                         # "20000"
# gtex_path                            <- args$gtex_path                           # "resources/GTEx/GTEx_different_versions_TCGA_cohort_map.RData"
# gtex_v8_path                         <- args$gtex_v8_path                        # "resources/GTEx/Tissue_specific_genes_based_on_GTEx_v8.RData"
# tissue_specific_proteins_path        <- args$tissue_specific_proteins_path       # "resources/GTEx/Tissue_specific_proteins_based_on_eGTEx_TSscore.RData"
# gtex_tissue_specific                 <- args$gtex_tissue_specific                # "resources/GTEx/GTEx_rna_protein_level_differentdatasets_per_bins_Tissue_specific.RData"
# gtex_output_path                     <- args$gtex_output_path                    # "resources/GTEx/GTEx_rna_protein_level_differentdatasets_per_bins.RData"

backbone_table                       <- "features_df.tsv"
window_size                          <- "20000"
gtex_path                            <- "resources/GTEx/GTEx_different_versions_TCGA_cohort_map.RData"
gtex_v8_path                         <- "resources/GTEx/Tissue_specific_genes_based_on_GTEx_v8.RData"
tissue_specific_proteins_path        <- "resources/GTEx/Tissue_specific_proteins_based_on_eGTEx_TSscore.RData"
gtex_tissue_specific                 <- "resources/GTEx/GTEx_rna_protein_level_differentdatasets_per_bins_Tissue_specific.RData"
gtex_output_path                     <- "resources/GTEx/GTEx_rna_protein_level_differentdatasets_per_bins.RData"

load(gtex_path)
common.tissues <- setdiff(colnames(datasets$eGTEX.rna),c("Gene.stable.ID","Gene.name"))

for(name in names(datasets)[-1]){
    common.tissues <- intersect(common.tissues, setdiff(colnames(datasets[[name]]),c("Gene.stable.ID","Gene.name")))}

gene.lists.per.bin <- function(bin){
  gene.set <- list("ppis.trans" = as.character(strsplit(backbone[backbone$binID == bin, "PPIs.trans.HIPPIE"],"::")[[1]]), #PPIs in trans
                   "all.int.trans" = union(as.character(strsplit(backbone[backbone$binID == bin, "Partners.trans"],"::")[[1]]),
                                           as.character(strsplit(backbone[backbone$binID == bin, "PPIs.trans.HIPPIE"],"::")[[1]])) #All interactions
  ) #For each bin, define a new gene sets
  names(gene.set) <- paste(bin, names(gene.set), sep = "_")
  return(gene.set)}

expression.scores <- function(genes){
  df <- exp.data[exp.data$Gene.name %in% genes,]
  if(nrow(df) == 0){ #In the case of no partner, PPIs or no covarage by expression data
    res <- data.frame("Cohorts" = colnames(df)[-c(1,2)], "Mean" = rep(NA, length(colnames(df)[-c(1,2)])), "Median" = rep(NA, length(colnames(df)[-c(1,2)])))
    }else{
    res <- data.frame("Cohorts" = colnames(df)[-c(1,2)], "Mean" = colMeans(df[,-c(1,2)], na.rm = TRUE), "Median" = apply(df[,-c(1,2)],2,median, na.rm = TRUE))
    }
  return(res)}

backbone <- read.table(file = backbone_table, header = T, sep = "\t")
backbone <- backbone[backbone$n_genes != 0,]
#selected.columns <- c("binID","chr","start_bin","end_bin", "genes","n_genes","Partners.trans","PPIs.trans")
#backbone <- backbone[, colnames(backbone) %in% selected.columns]

# non tissue specific

Results <- list()
gene.sets <- mclapply(backbone$binID, gene.lists.per.bin, mc.cores = 20)
gene.sets.merged <- list()

for(i in 1:length(gene.sets)){
    gene.sets.merged <- append(gene.sets.merged, gene.sets[[i]])
}

print(names(datasets))

for(dataset in names(datasets)){

    exp.data <- datasets[[dataset]]
    exp.data  <- exp.data[,colnames(exp.data) %in% c("Gene.stable.ID","Gene.name",common.tissues)]
    res.df <- mclapply(gene.sets.merged, expression.scores, mc.cores = 20)
    final.df <- do.call(rbind,res.df)

    if(window_size == "Arm"){
        final.df$binID <- unlist(lapply(rownames(final.df), function(x) strsplit(x,"_")[[1]][1]))
        final.df$GeneSet <- unlist(lapply(rownames(final.df), function(x) strsplit(x,"_")[[1]][2]))
    }else{
        final.df$binID <- unlist(lapply(rownames(final.df), function(x) paste(strsplit(x,"_")[[1]][1],strsplit(x,"_")[[1]][2],sep = "_")))
        final.df$GeneSet <- unlist(lapply(rownames(final.df), function(x) strsplit(x,"_")[[1]][3]))
    }

    final.df$GeneSet <- str_remove(final.df$GeneSet, paste0(".", paste(common.tissues,collapse = "|.")))
    final.df$GeneSet <- gsub("\\.\\d+","",final.df$GeneSet)
    final.df$GTExData <- dataset

    Results[[dataset]] <- final.df

}

save(Results, file = gtex_output_path)