options(scipen=999)                                                               # unable scientific notation

# Per-bin normal-tissue expression features (PENNE, "Features of normal tissue"): for each
# gene set of a bin, the median tissue-relative expression (z-score) of its genes in the
# GTEx cohort matched to the tumour type. NA when no gene of the set has expression data.
#   GTEx_genes, GTEx_partners_trans, GTEx_PPIs_trans, GTEx_all_int_trans
# Gene matching is by symbol (CORUM and Interactome INSIDER are symbol-based). GTEx symbols
# are first harmonised to the GENCODE release used for the backbone through the Ensembl gene
# ID; GTEx genes absent from that GENCODE annotation keep their GTEx symbol.

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_gene_sets   <- args$input_path_gene_sets
input_path_expression  <- args$input_path_expression
input_path_gencode     <- args$input_path_gencode
GTEx_cohort            <- args$GTEx_cohort
output_path            <- args$output_path

# quote = "" / comment.char = "": gene symbols can contain apostrophes (see compute_interactor_features.R)
gene_sets <- read.delim(input_path_gene_sets, stringsAsFactors = FALSE, quote = "", comment.char = "")
expression <- read.delim(input_path_expression, check.names = FALSE, stringsAsFactors = FALSE, quote = "", comment.char = "")
if(!GTEx_cohort %in% colnames(expression)){
  stop(paste0("GTEx cohort '", GTEx_cohort, "' not found in ", input_path_expression,
              " (available: ", paste(setdiff(colnames(expression), c("gene_id", "gene_name")), collapse = ", "), ")"))
}

gencode <- read.table(input_path_gencode, sep = "\t", header = F, colClasses = "character", quote = "", comment.char = "")[, 4:5]
colnames(gencode) <- c("gene_id", "gene_name")
gencode$gene_id <- sub("\\..*", "", gencode$gene_id)                               # strip version (and _PAR_Y)
gencode <- gencode[!duplicated(gencode$gene_id), ]

idx <- match(expression$gene_id, gencode$gene_id)
symbol <- ifelse(is.na(idx), expression$gene_name, gencode$gene_name[idx])
message(sprintf("GTEx genes: %d; matched to GENCODE by ID: %d; symbol changed by harmonisation: %d",
                nrow(expression), sum(!is.na(idx)), sum(!is.na(idx) & symbol != expression$gene_name)))

values <- data.frame(gene = symbol, value = expression[[GTEx_cohort]], stringsAsFactors = FALSE)
values <- values[!is.na(values$value), ]                                            # PENNE: median(na.rm = TRUE)

x <- merge(gene_sets, values, by = "gene")                                          # all rows sharing a symbol count (as in PENNE)

set_columns <- c(genes = "GTEx_genes", partners_trans = "GTEx_partners_trans",
                 ppis_trans = "GTEx_PPIs_trans", all_int_trans = "GTEx_all_int_trans")
bins <- unique(gene_sets$binID)
out <- data.frame(binID = bins, stringsAsFactors = FALSE)
for(set in names(set_columns)){
  sub <- x[x$gene_set == set, ]
  med <- if(nrow(sub) > 0) tapply(sub$value, sub$binID, median) else numeric(0)
  out[[set_columns[[set]]]] <- as.numeric(med[bins])                                # NA if no expressed gene
}

write.table(x = out, file = output_path, col.names = T, row.names = F, quote = F, sep = "\t")