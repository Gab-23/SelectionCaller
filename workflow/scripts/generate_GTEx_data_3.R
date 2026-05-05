# this shitty code is fortunately not mine
source("workflow/scripts/helpers/gene_lists_per_bin.R")
source("workflow/scripts/helpers/expression_scores.R")

options(scipen=999)                                                              # unable scientific notation     

library(parallel)

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)               # read args

backbone_table                       <- args$backbone_table                      # ""
window_size                          <- args$window_size
gtex_path                            <- args$gtex_path                           # "resources/GTEx/GTEx_different_versions_TCGA_cohort_map.RData"
gtex_v8_path                         <- args$gtex_v8_path                        # "resources/GTEx/Tissue_specific_genes_based_on_GTEx_v8.RData"
tissue_specific_proteins_path        <- args$tissue_specific_proteins_path       # "resources/GTEx/Tissue_specific_proteins_based_on_eGTEx_TSscore.RData"
gtex_tissue_specific                 <- args$gtex_tissue_specific                # "resources/GTEx/GTEx_rna_protein_level_differentdatasets_per_bins_Tissue_specific.RData"
gtex_output_path                     <- args$gtex_output_path                    # "resources/GTEx/GTEx_rna_protein_level_differentdatasets_per_bins.RData"

load(gtex_path)
common.tissues <- setdiff(colnames(datasets$eGTEX.rna),c("Gene.stable.ID","Gene.name"))
for(name in names(datasets)[-1]){
    common.tissues <- intersect(common.tissues, setdiff(colnames(datasets[[name]]),c("Gene.stable.ID","Gene.name")))}

backbone <- read.table(file = backbone_table, header = T, sep = "\t")
selected.columns <- c("bin","chr","start_bin","end_bin", "genes","n_genes","Partners.trans","PPIs.trans")

backbone <- backbone[, colnames(backbone) %in% selected.columns]
backbone <- backbone[backbone$n_genes != 0,]

# non tissue specific

Results <- list()
gene.sets <- mclapply(backbone$binID, gene.lists.per.bin, mc.cores = 20)
gene.sets.merged <- list()

for(i in 1:length(gene.sets)){
    gene.sets.merged <- append(gene.sets.merged, gene.sets[[i]])
}

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


# tissue-specific

rm(list = setdiff(ls(),c("backbone_tables_w_features_non_tissue_specific", "datasets", "common.tissues", "selected.columns")))
                                        
load(gtex_v8_path)                  #Tissue-specific genes defined based on RNA level data
load(tissue_specific_proteins_path) #Tissue-specific proteins defined based on eGTEx paper TS score

Results.TS <- list()

backbone <- read.table(file = backbone_table, header = T, sep = "\t")
backbone <- backbone[, colnames(backbone) %in% selected.columns]
backbone <- backbone[backbone$n_genes != 0,]

for(dataset in names(datasets)){

    if(dataset == "eGTEX.protein"){
        TS.list <- tissue.specific.proteins
    }else{ 
        TS.list <- tissue.specific.genes
    }

    tissue.final  <- c()
    for(tissue in common.tissues){
    
        tissue.genes <- TS.list[[tissue]]$Gene.name
        gene.sets <- mclapply(backbone$binID, gene.lists.per.bin)
        
        gene.sets.merged <- list()
        for(i in 1:length(gene.sets)){
            gene.sets.merged <- append(gene.sets.merged, gene.sets[[i]])
        }
        
        exp.data <- datasets[[dataset]]
        exp.data  <- exp.data[,colnames(exp.data) %in% c("Gene.stable.ID","Gene.name",tissue)]
        res.df <- mclapply(gene.sets.merged, expression.scores, mc.cores = 20)
        final.df <- do.call(rbind,res.df)
        
        if(window_size == "Arm"){
            final.df$bin <- unlist(lapply(rownames(final.df), function(x) strsplit(x,"_")[[1]][1]))
            final.df$GeneSet <- unlist(lapply(rownames(final.df), function(x) strsplit(x,"_")[[1]][2]))
        }else{
            final.df$bin <- unlist(lapply(rownames(final.df), function(x) paste(strsplit(x,"_")[[1]][1],strsplit(x,"_")[[1]][2],sep = "_")))
            final.df$GeneSet <- unlist(lapply(rownames(final.df), function(x) strsplit(x,"_")[[1]][3]))
        }

        final.df$GeneSet <- gsub("\\.\\d+","",final.df$GeneSet)
        final.df$GTExData <- dataset
        tissue.final <- rbind(tissue.final,final.df)

    }
    Results.TS[[dataset]] <- tissue.final
}

save(Results.TS, file = gtex_tissue_specific)