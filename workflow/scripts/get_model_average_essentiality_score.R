source("workflow/scripts/helpers/parse_depmap_colnames.R")

options(scipen=999)

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

input_path_scores       <- args$input_path_scores
input_path_genes        <- args$input_path_genes
output_path             <- args$output_path

crispr <- read.delim(input_path_scores, sep = ",")
crispr[is.na(crispr)] <- 0
colnames(crispr) <- parse_depmap_gene_symbols(colnames(crispr))

ess_genes <- read.table(input_path_genes, header = F, sep = "\t")
colnames(ess_genes) <- c("chr", "start", "end", "gene")

missing_genes <- setdiff(ess_genes$gene, colnames(crispr))
if (length(missing_genes) > 0) {
    warning(sprintf(
        "get_model_average_essentiality_score.R: %d gene(s) unexpectedly missing from scores columns: %s",
        length(missing_genes), paste(missing_genes, collapse = ", ")
    ))
}

crispr.f <- crispr[, c("ModelID", intersect(ess_genes$gene, colnames(crispr)))]

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