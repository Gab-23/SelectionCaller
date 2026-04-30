packages <- c("maditr")

installed <- rownames(installed.packages())
for (pkg in packages) {
  if (!pkg %in% installed) {
    install.packages(pkg, dependencies = TRUE)
  }
}

lapply(packages, library, character.only = TRUE)

GTEx_path                                   <- "./Data/GTEx_rna_protein_level_differentdatasets_per_bins.RData"
GTEx_ts_path                                <- "./Data/GTEx_rna_protein_level_differentdatasets_per_bins_Tissue_specific.RData"
backbone_path_with_features                 <- "./Data/Backbone_tables_with_all_features_GS_except_GTEx.RData"
backbone_outpath                            <- "./Data/Backbone_tables_with_all_features_GS.RData"
protein_levels_per_bins_path                <- "./Data/GTEx_rna_protein_level_differentdatasets_per_bins_HIPPIE.RData"
protein_levels_per_bins_ts_path             <- "./Data/GTEx_rna_protein_level_differentdatasets_per_bins_Tissue_specific_HIPPIE.RData"
protein_levels_cutoff_073_path              <- "./Data/GTEx_rna_protein_level_differentdatasets_per_bins_HIPPIE_cutoff_073.RData"
dataset <- "scaled.GTEx.v8.TPM"

load(GTEx_path) # Results
load(GTEx_ts_path) # Results.TS

# Dataset to be used
gtex.cohorts <- unique(Results$`1Mbp`$scaled.GTEx.v8.TPM$Cohorts)

# Feature data with the features previously added 
load(backbone_path_with_features)

levels <- names(backbone_tables_w_all_features_GS)

backbone_tables_w_all_features_GS_complete <- list()

for(level in levels){
  cohorts <- names(backbone_tables_w_all_features_GS[[level]])
  for(cohort in cohorts){
    
    m <- backbone_tables_w_all_features_GS[[level]][[cohort]]
    
    if(cohort %in% gtex.cohorts){
      
      # GTEx data (Some rows could be missing since for bins without any gene, there is no GTEx feature, assign NA to those bins)
      gtex <- Results[[level]][[dataset]]
      gtex <- gtex[gtex$Cohorts == cohort,]
      gtex <- gtex %>% reshape2::dcast(bin ~ GeneSet, value.var = "Median")
      colnames(gtex)[which(colnames(gtex) == "genes")] <- "genes.bin"
      colnames(gtex)[which(colnames(gtex) == "all.int.trans")] <- "all.int.trans_IntINSIDER"
      colnames(gtex)[which(colnames(gtex) == "ppis.trans")] <- "ppis.trans_IntINSIDER"
      # GTEx data (Some rows could be missing since for bins without any gene, there is no GTEx feature, assign NA to those bins)
      gtexTS <- Results.TS[[level]][[dataset]]
      gtexTS <- gtexTS[gtexTS$Cohorts == cohort,]
      gtexTS <- gtexTS %>% reshape2::dcast(bin ~ GeneSet, value.var = "Median")
      colnames(gtexTS)[which(colnames(gtexTS) == "genes")] <- "genes.bin"
      colnames(gtexTS)[2:5] <- paste(colnames(gtexTS)[2:5],"TS",sep = "_")
      colnames(gtexTS)[which(colnames(gtexTS) == "all.int.trans_TS")] <- "all.int.trans_TS_IntINSIDER"
      colnames(gtexTS)[which(colnames(gtexTS) == "ppis.trans_TS")] <- "ppis.trans_TS_IntINSIDER"
      
      m <- merge(m, gtex, by = "bin", all.x = T)
      m <- merge(m, gtexTS, by = "bin", all.x = T)
      
      backbone_tables_w_all_features_GS_complete[[level]][[cohort]] <- m
      
    }
    else{
      backbone_tables_w_all_features_GS_complete[[level]][[cohort]] <- m
    }
    
  
    }
}

save(backbone_tables_w_all_features_GS_complete, file = backbone_outpath)

rm(list = ls())
gc()