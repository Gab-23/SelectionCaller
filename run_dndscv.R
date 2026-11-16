## BURNER - dNdScv (Martincorena et al. 2017) sugli STESSI dati del pilota, per il confronto per gene.
## Solo sostituzioni singole (come il nostro target), stessi campioni COSMIC; dNdScv esclude da se' i campioni con piu' di
## 3000 mutazioni codificanti (max_coding_muts_per_sample = 3000, lo stesso criterio che usiamo noi).
## Uscita: la tabella sel_cv di dNdScv (p e q per gene: pmis_cv, ptrunc_cv, pallsubs_cv, ...) e la stima globale.
##
## Uso (ambiente R con il pacchetto dndscv):
##   Rscript run_dndscv.R --events resources/COSMIC/mutations/coding/f_breast_carcinoma.filtered.bed \
##     --refcds <RefCDS GRCh38 .rda> --cv hg19 --out_dir analyses/dndscv_BRCA
## --refcds: RefCDS per GRCh38 (le coordinate COSMIC del pilota sono GRCh38). --cv: covariate di dNdScv
## ("hg19" = matrice di default, indicizzata per NOME del gene; "none" = senza covariate, cioe' il modello
## "esposizione + sinonime del gene" senza regressione sulle covariate). Oppure --cv_file <TSV> (gene_name + covariate)
## con --cv_label <nome>: covariate proprie, per esempio il nostro background (make_dndscv_covariate.py).
## NOTE (16/11/2026): --sm <model> = substitution model of dNdScv: 192r_3w (default, trinucleotide context, as before),
## 12r_3w (12 substitution types, no context), 2r_3w (transitions/transversions only). With a model other than the
## default the output files get the suffix _sm<model> (e.g. sel_cv_hg19_hg38_epigenome_pcawg_sm12r_3w.tsv), so the
## default outputs are never overwritten. Used to measure how much the sequence context contributes to dNdScv.
args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)
suppressMessages(library(dndscv))
dir.create(args$out_dir, showWarnings = FALSE, recursive = TRUE)
cv <- if (is.null(args$cv)) "hg19" else args$cv
m <- read.delim(args$events, quote = "", stringsAsFactors = FALSE, colClasses = "character")
m <- m[m$GENOME_START == m$GENOME_STOP & toupper(m$GENOMIC_WT_ALLELE) %in% c("A","C","G","T") &
       toupper(m$GENOMIC_MUT_ALLELE) %in% c("A","C","G","T"), ]                               # sostituzioni singole
mut <- unique(data.frame(sampleID = m$COSMIC_SAMPLE_ID, chr = sub("^chr", "", m$CHROMOSOME),
                         pos = as.numeric(m$GENOME_START), ref = toupper(m$GENOMIC_WT_ALLELE),
                         mut = toupper(m$GENOMIC_MUT_ALLELE), stringsAsFactors = FALSE))         # una riga per evento
message(sprintf("sostituzioni singole uniche (campione x mutazione): %d, campioni: %d", nrow(mut), length(unique(mut$sampleID))))
## --cv_file: covariata/e da file (TSV: gene_name + una o piu' colonne). dNdScv vuole una riga per OGNI gene del suo RefCDS:
## i geni assenti dal file, o con valore mancante, ricevono la mediana della colonna (numero dichiarato nel log).
if (!is.null(args$cv_file)) {
    tab <- read.delim(args$cv_file, stringsAsFactors = FALSE)
    if (args$refcds == "hg38") data("refcds_GRCh38_hg38", package = "dndscv") else load(args$refcds)
    genes <- sapply(RefCDS, function(x) x$gene_name)
    M <- as.matrix(tab[match(genes, tab$gene_name), setdiff(colnames(tab), "gene_name"), drop = FALSE])
    rownames(M) <- genes
    for (j in seq_len(ncol(M))) {
        miss <- is.na(M[, j]); M[miss, j] <- median(M[, j], na.rm = TRUE)
        message(sprintf("covariata %s: %d geni del RefCDS, %d presenti nel file, %d riempiti con la mediana",
                        colnames(M)[j], length(genes), sum(!miss), sum(miss)))
    }
    cv_arg <- M; cv <- if (is.null(args$cv_label)) "custom" else args$cv_label
} else {
    cv_arg <- if (cv == "none") NULL else cv
}
## NOTE (16/11/2026): substitution model (--sm), checked against the models shipped with dNdScv
sm <- if (is.null(args$sm)) "192r_3w" else args$sm
allowed_sm <- c("192r_3w", "12r_3w", "2r_3w")
if (!(sm %in% allowed_sm)) stop(sprintf("--sm %s is not supported: use one of %s", sm, paste(allowed_sm, collapse = ", ")))
if (sm != "192r_3w") cv <- paste0(cv, "_sm", sm)          # suffix in the output names: default outputs untouched
message(sprintf("substitution model: %s; output label: %s", sm, cv))
res <- dndscv(mut, refdb = args$refcds, cv = cv_arg, sm = sm, max_coding_muts_per_sample = 3000)
write.table(res$sel_cv, file.path(args$out_dir, paste0("sel_cv_", cv, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(res$globaldnds, file.path(args$out_dir, paste0("globaldnds_", cv, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf("geni con q (tutte le sostituzioni) < 0.1: %d", sum(res$sel_cv$qallsubs_cv < 0.1, na.rm = TRUE)))
message(paste(head(res$sel_cv$gene_name[order(res$sel_cv$pallsubs_cv)], 20), collapse = ", "))