options(scipen=999)                                                               # unable scientific notation           

# args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

# input_path_backbone                                 <- args$input_path_backbone   # "resources/backbone/hg38_binned_20000_0.filtered.sorted.annotated.bed"
# ensembl_corum_path                                  <- args$ensembl_corum_path    # "resources/GTEx/Ensembl_CORUM_features_tissue_general.RData"
# hippie_path                                         <- args$hippie_path           # "resources/GTEx/HIPPIE_PPIs_tissue_general.RData"
# output_path                                         <- args$output_path           # "features_df.tsv"

input_path_backbone                                 <- "resources/backbone/hg38_binned_20000_0.filtered.sorted.annotated.bed"
ensembl_corum_path                                  <- "resources/GTEx/Ensembl_CORUM_features_tissue_general.RData"
hippie_path                                         <- "resources/GTEx/HIPPIE_PPIs_tissue_general.RData"
output_path                                         <- "features_df.tsv"

backbone_annotated <- read.table(input_path_backbone, sep = "\t", header = F)
colnames(backbone_annotated) <- c("chr", "start", "end", "binID", "genes")

load(ensembl_corum_path)
gene.table <- gene.table[,c(1,2,17:21)]

load(hippie_path)
colnames(gene.table.hippie) <- c("Gene.name","PPIs-HIPPIE","n_PPIs-HIPPIE", "PPIs_063-HIPPIE","n_PPIs_063-HIPPIE", "PPIs_073-HIPPIE","n_PPIs_073-HIPPIE")
gene.table.hippie[,c(3,5,7)] <- sapply(gene.table.hippie[,c(3,5,7)], as.numeric)

gene.table <- merge(gene.table, gene.table.hippie, by = "Gene.name", all.x = T)

backbone_tables_w_features <- list()

features.df <- c()

for(binID in backbone_annotated$binID){
    genes <- strsplit(backbone_annotated[backbone_annotated$binID == binID,"genes"],",")[[1]]
    df.genes <- gene.table[gene.table$Gene.name %in% genes,] # added this in case my annotation does not contain all genes from the other annotation. che schifo
    if(length(genes) == 0 || nrow(df.genes) == 0){
        row <- c(binID,"",length(genes),rep(NA,17))
    } else{

        partners <- unique(unlist(lapply(df.genes$Partners, function(x) strsplit(x,"::")[[1]])))
        partners <- partners[complete.cases(partners)]
        partners.trans <- setdiff(partners, genes)

        ppis <- unique(unlist(lapply(df.genes$`PPIs-HIPPIE`, function(x) strsplit(x,"::")[[1]])))
        ppis <- ppis[complete.cases(ppis)]
        ppis.trans <- setdiff(ppis,genes)

        ppis.063 <- unique(unlist(lapply(df.genes$`PPIs_063-HIPPIE`, function(x) strsplit(x,"::")[[1]])))
        ppis.063 <- ppis.063[complete.cases(ppis.063)]
        ppis.063.trans <- setdiff(ppis.063,genes)

        ppis.073 <- unique(unlist(lapply(df.genes$`PPIs_073-HIPPIE`, function(x) strsplit(x,"::")[[1]])))
        ppis.073 <- ppis.073[complete.cases(ppis.073)]
        ppis.073.trans <- setdiff(ppis.073,genes)

        row <- c(binID,
                paste(genes,collapse = "::"),
                length(genes),
                length(which(df.genes$ComplexSubunit == "True")),
                paste(partners, collapse = "::"),length(partners), 
                paste(partners.trans, collapse = "::"),length(partners.trans),
                paste(ppis, collapse = "::"),length(ppis),
                paste(ppis.trans, collapse = "::"), length(ppis.trans),
                paste(ppis.063, collapse = "::"),length(ppis.063),
                paste(ppis.063.trans, collapse = "::"), length(ppis.063.trans),
                paste(ppis.073, collapse = "::"),length(ppis.073),
                paste(ppis.073.trans, collapse = "::"), length(ppis.073.trans)
                )
        }

    features.df <- rbind(features.df,row)
    
    }

features.df <- as.data.frame(features.df)

colnames(features.df) <- c("binID", "genes", "n_genes", "total.complex.subunit",
                            "Partners","total_n_partners", "Partners.trans","total_n_partners.trans", "PPIs.HIPPIE","total_n_PPIs.HIPPIE",
                            "PPIs.trans.HIPPIE","total_n_PPIs.trans.HIPPIE", "PPIs_063.HIPPIE","total_n_PPIs_063.HIPPIE", 
                            "PPIs_063.trans.HIPPIE","total_n_PPIs_063.trans.HIPPIE", "PPIs_073.HIPPIE","total_n_PPIs_073.HIPPIE", 
                            "PPIs_073.trans.HIPPIE","total_n_PPIs_073.trans.HIPPIE")

features.df[,c(3,4,6,8,10,12,14,16,18,20)] <- sapply(features.df[,c(3,4,6,8,10,12,14,16,18,20)], as.numeric)

write.table(
    x = features.df,
    file = output_path,
    col.names = T,
    row.names = F,
    quote = F,
    sep = "\t"
)