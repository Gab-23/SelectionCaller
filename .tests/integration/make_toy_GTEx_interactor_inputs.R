# Builds toy inputs for the GTEx / protein-interaction module of the integration test,
# using the protein-coding genes of the toy GENCODE annotation already among the fixtures
# (so the toy backbone's genes actually have interactors and expression values):
#   resources/GTEx/v8/gene_median_tpm.gct.gz                    (all mapped GTEx tissues)
#   resources/protein_interactions/origin/toy_CORUM.RData       (object human_proteincomplex)
#   resources/protein_interactions/origin/toy_InteractomeINSIDER.RData (HS_interfaces_genenames)
# Usage (from repo root): Rscript .tests/integration/make_toy_GTEx_interactor_inputs.R \
#   <gencode gff3.gz> <GTEx->TCGA mapping> <output resources dir>

args <- commandArgs(trailingOnly = TRUE)
gff3_path <- args[1]; mapping_path <- args[2]; out_dir <- args[3]
set.seed(11)

lines <- readLines(gzfile(gff3_path))
lines <- lines[!startsWith(lines, "#")]
fields <- strsplit(lines, "\t", fixed = TRUE)
is_pc_gene <- vapply(fields, function(f) length(f) >= 9 && f[3] == "gene" &&
                       grepl("(^|;)gene_type=protein_coding(;|$)", f[9]), logical(1))
attr <- vapply(fields[is_pc_gene], `[`, character(1), 9)
get_attr <- function(a, key) sub(paste0(".*(^|;)", key, "=([^;]*).*"), "\\2", a)
genes <- unique(data.frame(id = get_attr(attr, "ID"), name = get_attr(attr, "gene_name"),
                           stringsAsFactors = FALSE))
# >= 2 genes is enough: one complex with all of them + a PPI between consecutive genes
if(nrow(genes) < 2) stop(paste("toy GENCODE has only", nrow(genes), "protein-coding genes; need >= 2"))
message(nrow(genes), " protein-coding toy genes")

# CORUM: one complex with all toy genes (every gene is a subunit with partners, whatever the
# number of toy genes), plus a single-subunit complex
human_proteincomplex <- data.frame(
  ComplexID = c(1, 2), ComplexName = c("toy complex", "toy monomer"),
  subunits.Gene.name. = c(paste(genes$name, collapse = ";"), genes$name[1]),
  stringsAsFactors = FALSE)
# Interactome INSIDER: a chain gene_i - gene_(i+1), plus a homodimer to be filtered out
n <- nrow(genes)
HS_interfaces_genenames <- data.frame(P1 = c(genes$name[-n], genes$name[1]),
                                      P2 = c(genes$name[-1], genes$name[1]),
                                      P1_IRES = "[]", P2_IRES = "[]", stringsAsFactors = FALSE)
dir.create(file.path(out_dir, "protein_interactions/origin"), recursive = TRUE, showWarnings = FALSE)
save(human_proteincomplex, file = file.path(out_dir, "protein_interactions/origin/toy_CORUM.RData"))
save(HS_interfaces_genenames, file = file.path(out_dir, "protein_interactions/origin/toy_InteractomeINSIDER.RData"))

# GTEx GCT v1.2: every tissue of the real mapping table (column names already in the
# make.names() form process_GTEx.R matches against), random positive TPMs
tissues <- unique(read.delim(mapping_path, stringsAsFactors = FALSE)$GTEx_tissue)
tpm <- matrix(round(rexp(n * length(tissues), rate = 0.1), 3), nrow = n, dimnames = list(NULL, tissues))
gct <- data.frame(Name = genes$id, Description = genes$name, tpm, check.names = FALSE)
dir.create(file.path(out_dir, "GTEx/v8"), recursive = TRUE, showWarnings = FALSE)
con <- gzfile(file.path(out_dir, "GTEx/v8/gene_median_tpm.gct.gz"), "w")
writeLines(c("#1.2", paste(n, length(tissues), sep = "\t")), con)
write.table(gct, con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)