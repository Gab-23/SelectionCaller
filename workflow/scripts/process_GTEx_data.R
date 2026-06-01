options(scipen=999)                                                              # unable scientific notation

library(dplyr)
library(reshape2)

# args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)               # read args
# 
# dataset                                     <- args$dataset                      # "scaled.GTEx.v8.TPM"
# GTEx_path                                   <- args$GTEx_path                    # "resources/GTEx/GTEx_rna_protein_level_differentdatasets_per_bins.RData"
# outname                                     <- args$outname                      # "prova_BRCA"
# tumor_type                                  <- args$tumor_type                   # "BRCA"
# outpath                                     <- args$outpath                      # paste0("data/GTEx/", dataset, "_", outname)

dataset                                     <- "scaled.GTEx.v8.TPM"
GTEx_path                                   <- "resources/GTEx/GTEx_rna_protein_level_differentdatasets_per_bins.RData"
outname                                     <- "prova_BRCA"
tumor_type                                  <- "BRCA"
outpath                                     <- paste0("data/GTEx/", dataset, "_", outname)

dataset <- "scaled.GTEx.v8.TPM"

load(GTEx_path) # Results

gtex.cohorts <- unique(Results[[dataset]]$Cohorts)

GTEx_list <- list()

for(cohort in gtex.cohorts){
  
  if(cohort %in% gtex.cohorts){
    
    gtex <- Results[[dataset]]
    gtex <- gtex[gtex$Cohorts == cohort,]
    gtex <- gtex %>% reshape2::dcast(binID ~ GeneSet, value.var = "Median")
    colnames(gtex)[which(colnames(gtex) == "genes")] <- "genes.bin"
    colnames(gtex)[which(colnames(gtex) == "all.int.trans")] <- "all.int.trans_IntINSIDER"
    colnames(gtex)[which(colnames(gtex) == "ppis.trans")] <- "ppis.trans_IntINSIDER"
    
    GTEx_list[[cohort]] <- gtex
    
  }
}

write.table(x = GTEx_list[[tissue_name]], 
            file = outpath, 
            header = T,
            row.names = F, 
            sep = "\t",
            quote = F)