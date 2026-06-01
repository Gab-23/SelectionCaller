options(scipen=999)     

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

input_path_scores       <- args$input_path_scores # "resources/DepMap/scores/BRCA_Female"
input_path_genes        <- args$input_path_genes #  "resources/DepMap/essential/BRCA_Female_0.5_4_v45.sorted.bed"
output_path             <- args$output_path # "resources/DepMap/score/BRCA_Female_0.5_4_v45.bed"

crispr <- read.delim(input_path_scores, sep = ",")
crispr[is.na(crispr)] <- 0
col_names <- gsub(pattern = "\\.\\..*", replacement = "", x = colnames(crispr))
col_names <- gsub(pattern = "\\.", replacement = "-", x = col_names)
col_names[1] <- "ModelID"
colnames(crispr) <- col_names

ess_genes <- read.table(input_path_genes, header = F, sep = "\t")
colnames(ess_genes) <- c("chr", "start", "end", "gene")

crispr.f <- crispr[, c("ModelID", ess_genes$gene[ess_genes$gene %in% colnames(crispr)])] # remove annotated fusion genes

mean_ess_score_across_models <- colMeans(crispr.f[,-1])

mean_ess_score_table <- data.frame(gene = names(mean_ess_score_across_models), 
								   ess_score = mean_ess_score_across_models)
								   
to.save <- merge(x = ess_genes, y = mean_ess_score_table, by = "gene")
to.save <- to.save[, c("chr", "start", "end", "gene", "ess_score")]

write.table(file = output_path, 
            x = to.save, 
            sep = "\t", 
            col.names = F,
            row.names = F, 
            quote = F)




