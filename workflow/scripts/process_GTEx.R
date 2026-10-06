options(scipen=999)                                                               # unable scientific notation

# Tissue-relative normal expression (PENNE, "Features of normal tissue"):
#   1. GTEx median gene TPM per tissue (GCT v1.2), Ensembl version suffix stripped;
#      duplicated IDs after stripping (PAR_Y copies) collapsed by maximum (as in PENNE).
#   2. GTEx tissues -> TCGA cohorts via the mapping table; cohorts with >1 tissue get the
#      median across tissues (as in PENNE's code).
#   3. Per gene, z-score across ALL mapped cohort columns (as in PENNE: a tissue mapped to
#      several cohorts - Kidney Cortex, Lung, Uterus - is counted once per cohort).
# Output: one row per gene, one z-score column per TCGA cohort. Downstream rules select
# the cohort they need; the z-score itself always uses every mapped cohort.

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_gct     <- args$input_path_gct
input_path_mapping <- args$input_path_mapping
output_path        <- args$output_path

# --- GTEx matrix -------------------------------------------------------------------------
# skip=2: GCT version line + dimensions line. check.names=TRUE (read.delim default) turns
# "Breast - Mammary Tissue" into "Breast...Mammary.Tissue", the form used in the mapping
# table (PENNE Supplementary Table 9).
gct <- read.delim(gzfile(input_path_gct), skip = 2, check.names = TRUE, stringsAsFactors = FALSE,
                  quote = "", comment.char = "")                                  # symbols may contain apostrophes
if(!all(c("Name", "Description") %in% colnames(gct))){
  stop("GCT file must have 'Name' and 'Description' as first columns")
}
gct$Name <- sub("\\..*", "", gct$Name)                                           # strip Ensembl version (and _PAR_Y)

mapping <- read.delim(input_path_mapping, stringsAsFactors = FALSE, quote = "", comment.char = "")
if(!all(c("GTEx_tissue", "TCGA_cohort") %in% colnames(mapping))){
  stop("mapping file must have columns 'GTEx_tissue' and 'TCGA_cohort'")
}
mapping <- mapping[!is.na(mapping$TCGA_cohort) & mapping$TCGA_cohort != "", ]

# Fail loudly instead of silently dropping tissues (PENNE's merge() dropped them): a
# renamed tissue in a new GTEx release would otherwise change a cohort's definition.
missing_tissues <- setdiff(unique(mapping$GTEx_tissue), colnames(gct))
if(length(missing_tissues) > 0){
  stop(paste0("mapped GTEx tissues not found in the GCT file: ",
              paste(missing_tissues, collapse = ", ")))
}

tissues <- unique(mapping$GTEx_tissue)
expr <- as.matrix(gct[, tissues, drop = FALSE])
storage.mode(expr) <- "numeric"

# Collapse duplicated gene IDs by maximum (PENNE: group_by(Name) + max)
gene_ids <- unique(gct$Name)
if(any(duplicated(gct$Name))){
  expr <- do.call(rbind, lapply(split(seq_len(nrow(expr)), factor(gct$Name, levels = gene_ids)),
                                function(idx) apply(expr[idx, , drop = FALSE], 2, max, na.rm = TRUE)))
}
gene_names <- gct$Description[match(gene_ids, gct$Name)]                         # first symbol per ID

# --- tissues -> cohorts (median) ----------------------------------------------------------
cohorts <- sort(unique(mapping$TCGA_cohort))
cohort_expr <- sapply(cohorts, function(cohort){
  cohort_tissues <- unique(mapping$GTEx_tissue[mapping$TCGA_cohort == cohort])
  apply(expr[, cohort_tissues, drop = FALSE], 1, median, na.rm = TRUE)
})
cohort_expr <- matrix(cohort_expr, nrow = length(gene_ids), dimnames = list(gene_ids, cohorts))

# --- tissue-relative z-score per gene across cohorts ------------------------------------
# Same as scale() on each row (sample SD). Genes with SD = 0 (e.g. not expressed anywhere)
# get NA, as scale() would give NaN.
row_mean <- rowMeans(cohort_expr, na.rm = TRUE)
row_sd   <- apply(cohort_expr, 1, sd, na.rm = TRUE)
zscores  <- (cohort_expr - row_mean) / row_sd
zscores[!is.finite(zscores)] <- NA

out <- data.frame(gene_id = gene_ids, gene_name = gene_names, zscores,
                  check.names = FALSE, stringsAsFactors = FALSE)

write.table(
    x = out,
    file = output_path,
    col.names = T,
    row.names = F,
    quote = F,
    sep = "\t"
)