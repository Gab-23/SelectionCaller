# this shitty code is fortunately not mine

options(scipen=999)                                                              # unable scientific notation

library(dplyr)
library(reshape2)

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)               # read args

GTEx_path                                   <- args$GTEx_path                    # "resources/GTEx/GTEx_rna_protein_level_differentdatasets_per_bins.RData"
outpath                                     <- args$outpath
dataset                                     <- "scaled.GTEx.v8.TPM"

load(GTEx_path) # Results

gtex.cohorts <- unique(Results$scaled.GTEx.v8.TPM$Cohorts)

backbone_tables_w_all_features_GS_complete <- list()
cohorts <- names(backbone_tables_w_all_features_GS)

for(cohort in gtex.cohorts){
  
  m <- backbone_tables_w_all_features_GS[[cohort]]
  
  if(cohort %in% gtex.cohorts){
    
    gtex <- Results[[dataset]]
    gtex <- gtex[gtex$Cohorts == cohort,]
    gtex <- gtex %>% reshape2::dcast(bin ~ GeneSet, value.var = "Median")
    colnames(gtex)[which(colnames(gtex) == "genes")] <- "genes.bin"
    colnames(gtex)[which(colnames(gtex) == "all.int.trans")] <- "all.int.trans_IntINSIDER"
    colnames(gtex)[which(colnames(gtex) == "ppis.trans")] <- "ppis.trans_IntINSIDER"
    
    m <- merge(m, gtex, by = "binID", all.x = T)
    
    backbone_tables_w_all_features_GS_complete[[cohort]] <- m
    
  } else{
    backbone_tables_w_all_features_GS_complete[[cohort]] <- m
  }
}

save(backbone_tables_w_all_features_GS_complete, file = backbone_outpath)