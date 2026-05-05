# this shitty code is fortunately not mine

options(scipen=999)                                                                                # unable scientific notation

library(openxlsx)

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                                 # read args

info_path                                             <- args$info_path                            # "resources/GTEx/GTEX_tissuenames_tcganames.xlsx"
gtex_protein_path                                     <- args$gtex_protein_path                    # "resources/GTEx/Table_S1_gene_info_at_protein_level.xlsx"
gtex_mrna_path                                        <- args$gtex_mrna_path                       # "resources/GTEx/Table_S2_gene_info_at_RNA_level.xlsx"
tissuenames_tcganames_path                            <- args$tissuenames_tcganames_path           # "resources/GTEx/eGTEX_tissuenames_tcganames.xlsx"
gtex_protein_tissue_median_path                       <- args$gtex_protein_tissue_median_path      # "resources/GTEx/eGTEX_mrna_protein_tissue_median.RData"
mrna_v8_path                                          <- args$mrna_v8_path                         # "resources/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct"
gtex_outpath                                          <- args$gtex_outpath                         # "resources/GTEx/GTEx_different_versions_TCGA_cohort_map.RData"
tissue_specific_genes_v8                              <- args$tissue_specific_genes_v8             # "resources/GTEx/Tissue_specific_genes_based_on_GTEx_v8.RData"
tissuespec_proteins                                   <- args$tissuespec_proteins                  # "resources/GTEx/Tissue_specific_proteins_based_on_eGTEx_TSscore.RData"

info <- read.xlsx(info_path,startRow = 2)

gtex.protein <- read.xlsx(gtex_protein_path,sheet = 9)
gtex.mrna <- read.xlsx(gtex_mrna_path,sheet = 3)

geneID_name <- read.xlsx(gtex_protein_path,sheet = 10)
gtex.protein <- merge(geneID_name[,c(1,4)], gtex.protein, by.x="ensembl_id",by.y = "gene.id")
gtex.mrna <- merge(geneID_name[,c(1,4)], gtex.mrna, by.x="ensembl_id",by.y = "gene.id")

tissue.names <- read.xlsx(gtex_protein_path,sheet = 12, startRow = 2)
tissue.names <- rbind(tissue.names[,1:2],tissue.names[3:4])
tissue.names$Colnames <- colnames(gtex.protein)[3:34]

write.xlsx(tissue.names, file = tissuenames_tcganames_path)
save(gtex.mrna, gtex.protein, file = gtex_protein_tissue_median_path)

rm(list = setdiff(ls(),c("info","datasets")))

datasets <- list()

load(gtex_protein_tissue_median_path)

# eGTEX.rna
rownames(gtex.mrna) <- gtex.mrna$ensembl_id
m <- as.data.frame(t(gtex.mrna[,-c(1,2)]))
m$Colnames <- rownames(m)
m <- merge(info[,5:6],m,by="Colnames")
m <- m[complete.cases(m$TCGA),]
m <- m %>% group_by(TCGA) %>% summarise(across(where(is.numeric), ~ median(.x, na.rm = TRUE)))
m <- as.data.frame(m)
rownames(m) <- m$TCGA
m <- as.data.frame(t(m[,-1]))
m$Gene.stable.ID <- rownames(m)
m <- merge(gtex.mrna[,1:2],m,by.x = "ensembl_id" ,by.y = "Gene.stable.ID")
colnames(m)[1:2] <- c("Gene.stable.ID","Gene.name") 
datasets[["eGTEX.rna"]] <- m

#eGTEX protein
rownames(gtex.protein) <- gtex.protein$ensembl_id
m <- as.data.frame(t(gtex.protein[,-c(1,2)]))
m$Colnames <- rownames(m)
m <- merge(info[,5:6],m,by="Colnames")
m <- m[complete.cases(m$TCGA),]
m <- m %>% group_by(TCGA) %>% summarise(across(where(is.numeric), ~ median(.x, na.rm = TRUE)))
m <- as.data.frame(m)
rownames(m) <- m$TCGA
m <- as.data.frame(t(m[,-1]))
m$Gene.stable.ID <- rownames(m)
m <- merge(gtex.protein[,1:2],m,by.x = "ensembl_id" ,by.y = "Gene.stable.ID")
colnames(m)[1:2] <- c("Gene.stable.ID","Gene.name") 
datasets[["eGTEX.protein"]] <- m

rm(list = setdiff(ls(),c("info","datasets")))

#GTEx.v6.RPKM
mrna.v6p <- read.delim(gtex_v6_RNA_median,skip = 2)
mrna.v6p$Name <- unlist(lapply(mrna.v6p$Name, function(x) gsub("\\..*","",x)))
rownames(mrna.v6p) <- mrna.v6p$Name
m <- as.data.frame(t(mrna.v6p[,-c(1,2)]))
m$V6p <- rownames(m)
m <- merge(info[,c(1,6)],m,by="V6p")
m <- m[complete.cases(m$TCGA),]
m <- m %>% group_by(TCGA) %>% summarise(across(where(is.numeric), ~ median(.x, na.rm = TRUE)))
m <- as.data.frame(m)
rownames(m) <- m$TCGA
m <- as.data.frame(t(m[,-1]))
m$Gene.stable.ID <- rownames(m)
m <- merge(mrna.v6p[,1:2],m,by.x = "Name" ,by.y = "Gene.stable.ID")
colnames(m)[1:2] <- c("Gene.stable.ID","Gene.name")
datasets[["GTEx.v6.RPKM"]] <- m


#GTEx.v8.TPM
mrna.v8 <- read.delim(mrna_v8_path,skip = 2)
mrna.v8$Name <- unlist(lapply(mrna.v8$Name, function(x) gsub("\\..*","",x)))
gene.info <- mrna.v8[,1:2]
gene.info <- gene.info[!duplicated(gene.info$Name),]
mrna.v8 <- mrna.v8 %>% group_by(Name) %>% summarise(across(where(is.numeric), ~ max(.x, na.rm = TRUE))) #For duplicated ensembl IDs, keep the maximum
mrna.v8 <- merge(gene.info, mrna.v8, by = "Name")
rownames(mrna.v8) <- mrna.v8$Name
m <- as.data.frame(t(mrna.v8[,-c(1,2)]))
m$V8 <- rownames(m)
m <- merge(info[,c(1,6)],m,by.x="V6p",by.y = "V8")
m <- m[complete.cases(m$TCGA),]
m <- m %>% group_by(TCGA) %>% summarise(across(where(is.numeric), ~ median(.x, na.rm = TRUE)))
m <- as.data.frame(m)
rownames(m) <- m$TCGA
m <- as.data.frame(t(m[,-1]))
m$Gene.stable.ID <- rownames(m)
m <- merge(mrna.v8[,1:2],m,by.x = "Name" ,by.y = "Gene.stable.ID")
colnames(m)[1:2] <- c("Gene.stable.ID","Gene.name")
datasets[["GTEx.v8.TPM"]] <- m

save(datasets, file = gtex_outpath)

#scaled.GTEx.v8.TPM
m <- datasets$GTEx.v8.TPM
m.scaled <- c()
for(i in 1:nrow(m)){
  v <- as.numeric(m[i,-c(1,2)])
  v.scaled <- as.numeric(scale(v))
  m.scaled <- rbind(m.scaled, c(m$Gene.stable.ID[i],m$Gene.name[i],v.scaled))
  print(i)}
m.scaled <- as.data.frame(m.scaled)
colnames(m.scaled) <- colnames(m)
m.scaled[,3:25] <- sapply(m.scaled[,3:25], as.numeric)
datasets[["scaled.GTEx.v8.TPM"]] <- m.scaled

save(datasets, file = gtex_outpath)

#Tissue-specific genes
#First remove genes less than 1 TPM
#Then filter based on z-score (abs z-score > 2, tissue specific) - using absolute value is different from Uri-Ben David paper (preferential low expression could play a role in selection)
#I did not use eGTEx paper's TS because the number of the genes in RNA data is limited with proteome covarage
load(gtex_outpath)
m <- datasets$GTEx.v8.TPM
zscores <- datasets$scaled.GTEx.v8.TPM
tissue.specific.genes <- list()
cohorts <- colnames(m)[3:25]
for(cohort in cohorts){
  df <- m[,c("Gene.stable.ID","Gene.name",cohort)]
  colnames(df)[3] <- "tpm"
  df <- df[df$tpm > 1,]
  df.scaled <- zscores[,c("Gene.stable.ID","Gene.name",cohort)]
  colnames(df.scaled)[3] <- "zscore"
  df.scaled <- df.scaled[abs(df.scaled$zscore) > 2,]
  genes <- intersect(df$Gene.stable.ID,df.scaled$Gene.stable.ID)
  tissue.specific.genes[[cohort]] <- df[df$Gene.stable.ID %in% genes, 1:2]}

save(tissue.specific.genes, file = tissue_specific_genes_v8)

#Tissue-specific proteins
#Based on TS score in https://www.sciencedirect.com/science/article/pii/S0092867420310783?via%3Dihub - eGTEx paper
gtex.protein.TS <- read.xlsx(gtex_protein_path,sheet = 11) #12627 - Median relative protein abundance in each tissue type
for(col in 5:ncol(gtex.protein.TS)){
  new.values <- as.numeric(unlist(lapply(gtex.protein.TS[,col], function(x) strsplit(x,";")[[1]][1])))
  gtex.protein.TS[,col] <- new.values}

tissues <- colnames(gtex.protein.TS)[5:ncol(gtex.protein.TS)]

tissue.specific.proteins.tissues <- list()

for(tissue in tissues){
  df <- gtex.protein.TS[,c("ensembl_id","hgnc_symbol",tissue)]
  df <- df[complete.cases(df),]
  colnames(df)[3] <- "tissue"
  df <- df[abs(df$tissue) >= 2.5,]
  df <- df[,1:2]
  colnames(df) <- c("Gene.stable.ID","Gene.name")
  tissue.specific.proteins.tissues[[tissue]] <- df}

#Converting tissue names to TCGA cohorts
info <- read.xlsx(info_path,startRow = 2)
info <- info[,5:6]
info <- info[complete.cases(info),]
tissue.specific.proteins <- list()
for(cohort in unique(info$TCGA)){
  names <- info[info$TCGA == cohort,]$Colnames
  genes <- c()
  for(name in names){
    df <- tissue.specific.proteins.tissues[[name]]
    genes <- rbind(genes,df)}
  genes <- genes[!duplicated(genes$Gene.stable.ID),]
  tissue.specific.proteins[[cohort]] <- genes
}

save(tissue.specific.proteins, file = tissuespec_proteins)
