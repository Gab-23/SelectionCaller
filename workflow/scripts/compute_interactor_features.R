options(scipen=999)                                                               # unable scientific notation

# Per-bin interactor features (PENNE, "Features inferred from the gene table") and the
# per-bin gene sets reused by the GTEx features.
#   complex_density           = n. bin genes that are CORUM subunits / n. bin genes
#   n_complex_partners_trans  = n. unique CORUM partners of bin genes located outside the bin
#   n_PPIs_trans              = n. unique Interactome INSIDER partners located outside the bin
# "In trans" = not among the genes of the same bin. Bins without protein-coding genes are
# absent from the annotated backbone, so they get NA after the merge in combine_dataset;
# bins with genes but no partners get 0 (as in PENNE).
# Gene sets (long format: binID, gene_set, gene): genes, partners_trans, ppis_trans,
# all_int_trans (= union of the two trans sets).

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_backbone     <- args$input_path_backbone
input_path_interactors  <- args$input_path_interactors
output_path_features    <- args$output_path_features
output_path_gene_sets   <- args$output_path_gene_sets

# quote = "" and comment.char = "": gene symbols can contain apostrophes (e.g. Interactome
# INSIDER "PSG11s'"), which read.table would otherwise treat as quotes, silently merging rows.
backbone <- read.table(input_path_backbone, sep = "\t", header = F, colClasses = "character", quote = "", comment.char = "",
                       col.names = c("chr", "start", "end", "binID", "genes"))
genes_split <- strsplit(backbone$genes, ",")
bin_genes <- unique(data.frame(binID = rep(backbone$binID, lengths(genes_split)),
                               gene = unlist(genes_split), stringsAsFactors = FALSE))
bin_genes <- bin_genes[!is.na(bin_genes$gene) & bin_genes$gene != "", ]

interactors <- read.delim(input_path_interactors, stringsAsFactors = FALSE, quote = "", comment.char = "")

trans_partners <- function(src){
  pairs <- interactors[interactors$source == src & !is.na(interactors$partner), c("gene", "partner")]
  x <- merge(bin_genes, pairs, by = "gene")
  x <- unique(x[, c("binID", "partner")])
  colnames(x) <- c("binID", "gene")
  x[!paste(x$binID, x$gene) %in% paste(bin_genes$binID, bin_genes$gene), ]      # keep trans only
}

partners_trans <- trans_partners("CORUM")
ppis_trans     <- trans_partners("InteractomeINSIDER")
all_int_trans  <- unique(rbind(partners_trans, ppis_trans))

gene_sets <- rbind(data.frame(binID = bin_genes$binID,      gene_set = "genes",          gene = bin_genes$gene),
                   data.frame(binID = partners_trans$binID, gene_set = "partners_trans", gene = partners_trans$gene),
                   data.frame(binID = ppis_trans$binID,     gene_set = "ppis_trans",     gene = ppis_trans$gene),
                   data.frame(binID = all_int_trans$binID,  gene_set = "all_int_trans",  gene = all_int_trans$gene))

bins <- unique(bin_genes$binID)
count_per_bin <- function(df){ n <- table(factor(df$binID, levels = bins)); as.numeric(n[bins]) }

corum_subunits <- unique(interactors$gene[interactors$source == "CORUM"])
n_genes    <- count_per_bin(bin_genes)
n_subunits <- count_per_bin(bin_genes[bin_genes$gene %in% corum_subunits, ])

features <- data.frame(binID = bins,
                       complex_density = n_subunits / n_genes,
                       n_complex_partners_trans = count_per_bin(partners_trans),
                       n_PPIs_trans = count_per_bin(ppis_trans),
                       stringsAsFactors = FALSE)

write.table(x = features,  file = output_path_features,  col.names = T, row.names = F, quote = F, sep = "\t")
write.table(x = gene_sets, file = output_path_gene_sets, col.names = T, row.names = F, quote = F, sep = "\t")