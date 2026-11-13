options(scipen=999)                                                               # unable scientific notation

# Conteggi di mutazioni codificanti per bin, a livello di EVENTO (campione x mutazione).
# Diverso da compute_coding_mutations_overlap.R, che lavora sul file collassato
# (sort_collapse_mutations: `sort | uniq` sulle coordinate) e quindi conta ogni mutazione
# UNA volta sola, qualunque sia il numero di campioni che la portano.
#
# Colonne in output (una riga per bin con almeno una mutazione; gli altri bin diventano 0 in
# combine_dataset):
#   mut_n_events   eventi = coppie (campione, mutazione) distinte -> include la ricorrenza
#   mut_n_samples  campioni distinti con almeno una mutazione nel bin
#   mut_n_sites    mutazioni distinte (GENOMIC_MUTATION_ID) -> ricorrenza collassata
#   --- AGGIUNTE 22/10/2026: eventi per CLASSE, solo SOSTITUZIONI SINGOLE (convenzioni di dNdScv) ---
#   mut_n_syn            sinonime
#   mut_n_mis            missense
#   mut_n_trunc          troncanti: nonsense e siti di splicing essenziali (acceptor / donor)
#   mut_n_other_coding   altre codificanti (start_lost, stop_lost, stop_retained_variant, ...): tenute a
#                        parte, non assegnate a una classe (decisione rimandata)
# Le prime tre colonne sono calcolate ESATTAMENTE come prima (tutti gli eventi, qualsiasi conseguenza).
#
# Sostituzione singola: GENOME_START == GENOME_STOP e GENOMIC_WT_ALLELE / GENOMIC_MUT_ALLELE = una base ACGT.
# Classe di una mutazione (il file COSMIC ha una riga per trascritto), parametro transcript_mode:
#   canonical  conseguenza nei trascritti Ensembl_canonical (tag del GTF); se la mutazione tocca i canonici
#              di piu' geni (geni sovrapposti) vale la classe piu' grave; le mutazioni che non cadono in
#              nessun trascritto canonico non entrano nelle colonne per classe (numero riportato nel log)
#   all        classe piu' grave tra tutti i trascritti (una mutazione codificante in almeno un trascritto
#              conta come codificante: unione a livello di gene)
#   gravita':  troncante > missense > altra codificante > sinonima > non codificante
# Ogni evento conta UNA volta nel bin, con la sua classe.
#
# Campioni ipermutati: se hypermutator_max != "none", si escludono i campioni con piu' di hypermutator_max
# eventi CODIFICANTI (sostituzioni e indel; conseguenza codificante in almeno un trascritto), come in
# Martincorena et al. 2017 (> 3000 mutazioni codificanti). La soglia si applica PRIMA di tutti i conteggi.
# Una mutazione e' assegnata al bin che contiene GENOME_START (1-based).
#
# FIX 28/10/2026 - RIANNOTAZIONE degli eventi (parametro event_annotation, default "reference"):
#   PROBLEMA: con "cosmic" la classe viene dalle righe del file COSMIC sui trascritti della modalita'; ma COSMIC usa un
#   proprio insieme di trascritti, e per alcuni geni il canonico Ensembl v45 non c'e' MAI (es. AKT1, CASP8, RUNX1): tutti i
#   loro eventi sparivano dalle colonne per classe (8.7% degli eventi codificanti nella coorte BRCA, verificato).
#   "reference": ogni SOSTITUZIONE SINGOLA viene classificata da noi, sulla sequenza di riferimento e sui trascritti della
#   modalita', con le STESSE regole e lo STESSO codice di compute_coding_sites.R (helpers/coding_sites_lib.R): la classe
#   di un evento e' quella del sito (posizione + base mutata) nella tabella dei siti; eventi ed esposizione sono definiti
#   nello stesso modo per costruzione (come dNdScv, che annota da se' le mutazioni sul suo RefCDS).
#   - evento che non e' un sito codificante dei trascritti della modalita' (introne, UTR, fuori dai geni): nessuna classe;
#   - base di riferimento COSMIC diversa dal genoma: nessuna classe (numero nel log);
#   - indel e sostituzioni multi-base: come prima, nessuna classe;
#   - colonne mut_n_events / mut_n_samples / mut_n_sites e soglia degli ipermutati: INVARIATE (usano ancora le
#     conseguenze COSMIC su qualsiasi trascritto).
#   "cosmic": comportamento precedente, invariato.
# NOTE (28/10/2026): riepiloghi in R base (niente dplyr), per girare nello stesso ambiente di compute_coding_sites.R.

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_mutations    <- args$input_path_mutations       # COSMIC GenomeScreensMutant filtrato (con header)
input_path_backbone     <- args$input_path_backbone
input_path_gtf          <- args$input_path_gtf             # GENCODE (gz): trascritti con tag Ensembl_canonical
hypermutator_max        <- as.character(args$hypermutator_max)
transcript_mode         <- if (is.null(args$transcript_mode)) "canonical" else as.character(args$transcript_mode)
sample_col              <- if (is.null(args$sample_col)) "COSMIC_SAMPLE_ID" else as.character(args$sample_col)
output_path             <- args$output_path
output_path_samples     <- args$output_path_samples
event_annotation        <- if (is.null(args$event_annotation)) "reference" else as.character(args$event_annotation)
input_path_genome       <- args$input_path_genome          # FASTA (solo per event_annotation = "reference")
if (!event_annotation %in% c("reference", "cosmic")) stop("event_annotation deve essere 'reference' o 'cosmic', trovato: ", event_annotation)
if (event_annotation == "reference" && is.null(input_path_genome)) stop("event_annotation = 'reference' richiede --input_path_genome")

if (!transcript_mode %in% c("canonical", "all")) stop("transcript_mode deve essere 'canonical' o 'all', trovato: ", transcript_mode)

strip_chr <- function(x) sub("^chr", "", x)
strip_version <- function(x) sub("\\..*$", "", x)

# ---------------- classificazione delle conseguenze ----------------
TRUNC_TERMS  <- c("stop_gained", "splice_acceptor_variant", "splice_donor_variant")
MIS_TERMS    <- c("missense_variant")
OTHER_TERMS  <- c("start_lost", "stop_lost", "stop_retained_variant", "coding_sequence_variant",
                  "protein_altering_variant", "incomplete_terminal_codon_variant", "start_retained_variant")
SYN_TERMS    <- c("synonymous_variant")
INDEL_TERMS  <- c("frameshift_variant", "inframe_insertion", "inframe_deletion")   # codificanti, solo per gli ipermutati
CLASS_NAMES  <- c("noncoding", "syn", "other_coding", "mis", "trunc")             # indice = gravita' (1 = meno grave)

severity_of <- function(desc) {
    # gravita' di una MUTATION_DESCRIPTION (termini separati da virgola): 5 trunc, 4 mis, 3 other, 2 syn, 1 non codificante
    u <- unique(desc); sev <- integer(length(u))
    for (i in seq_along(u)) {
        t <- strsplit(u[i], ",", fixed = TRUE)[[1]]
        sev[i] <- if (any(t %in% TRUNC_TERMS)) 5L else if (any(t %in% MIS_TERMS)) 4L else
                  if (any(t %in% OTHER_TERMS)) 3L else if (any(t %in% SYN_TERMS)) 2L else 1L
    }
    sev[match(desc, u)]
}
is_coding_any <- function(desc) {
    # codificante in senso lato (sostituzioni e indel): per la soglia degli ipermutati
    vapply(strsplit(desc, ",", fixed = TRUE), function(t) any(t %in% c(TRUNC_TERMS, MIS_TERMS, OTHER_TERMS, SYN_TERMS, INDEL_TERMS)), logical(1))
}

# ---------------- mutazioni ----------------
muts <- read.delim(file = input_path_mutations, header = TRUE, sep = "\t", quote = "",
                   comment.char = "", colClasses = "character", check.names = FALSE)
required <- c(sample_col, "GENOMIC_MUTATION_ID", "CHROMOSOME", "GENOME_START", "GENOME_STOP",
              "TRANSCRIPT_ACCESSION", "MUTATION_DESCRIPTION", "GENOMIC_WT_ALLELE", "GENOMIC_MUT_ALLELE")
missing <- setdiff(required, colnames(muts))
if (length(missing) > 0) stop("colonne mancanti nel file COSMIC: ", paste(missing, collapse = ", "))

rows <- data.frame(sample = muts[[sample_col]],
                   mut_id = ifelse(muts$GENOMIC_MUTATION_ID == "" | is.na(muts$GENOMIC_MUTATION_ID),
                                   paste(muts$CHROMOSOME, muts$GENOME_START, muts$GENOME_STOP, sep = ":"),
                                   muts$GENOMIC_MUTATION_ID),
                   chr = strip_chr(muts$CHROMOSOME),
                   pos = suppressWarnings(as.numeric(muts$GENOME_START)),
                   stop = suppressWarnings(as.numeric(muts$GENOME_STOP)),
                   tx = strip_version(muts$TRANSCRIPT_ACCESSION),
                   desc = muts$MUTATION_DESCRIPTION,
                   wt = toupper(muts$GENOMIC_WT_ALLELE), mt = toupper(muts$GENOMIC_MUT_ALLELE),
                   stringsAsFactors = FALSE)
rm(muts)
rows <- rows[!is.na(rows$pos) & rows$sample != "", ]
rows$sev <- severity_of(rows$desc)
rows$snv <- !is.na(rows$stop) & rows$pos == rows$stop & rows$wt %in% c("A", "C", "G", "T") & rows$mt %in% c("A", "C", "G", "T")

# trascritti canonici dal GTF
if (event_annotation == "cosmic" && transcript_mode == "canonical") {
    gl <- readLines(gzfile(input_path_gtf))
    gl <- gl[grepl("\ttranscript\t", gl, fixed = TRUE) & grepl('tag "Ensembl_canonical"', gl, fixed = TRUE)]
    canonical <- strip_version(sub('.*transcript_id "([^"]+)".*', "\\1", gl))
    rm(gl)
    message(sprintf("trascritti Ensembl_canonical nel GTF: %d", length(canonical)))
    rows$in_mode <- rows$tx %in% canonical
} else {
    rows$in_mode <- TRUE
}

# un evento = (campione, mutazione); attributi dell'evento aggregati sulle sue righe (trascritti)
key <- paste(rows$sample, rows$mut_id, sep = "\r")
coding_any <- tapply(is_coding_any(rows$desc), key, any)                        # per la soglia degli ipermutati
sev_mode <- tapply(ifelse(rows$in_mode, rows$sev, 0L), key, max)                # 0 = nessun trascritto della modalita'
first <- !duplicated(key)
events <- data.frame(sample = rows$sample[first], mut_id = rows$mut_id[first], chr = rows$chr[first],
                     pos = rows$pos[first], snv = rows$snv[first], wt = rows$wt[first], mt = rows$mt[first], stringsAsFactors = FALSE)
k1 <- key[first]
events$coding_any <- as.logical(coding_any[k1]); events$sev <- as.integer(sev_mode[k1])
rm(rows)

# ---------------- ipermutati (eventi codificanti per campione) ----------------
per_sample <- data.frame(sample = sort(unique(events$sample)), stringsAsFactors = FALSE)
per_sample$n_events <- as.integer(table(factor(events$sample, levels = per_sample$sample)))
per_sample$n_coding_events <- as.integer(tapply(events$coding_any, factor(events$sample, levels = per_sample$sample), sum))
if (hypermutator_max %in% c("none", "None", "NULL", "")) {
    per_sample$kept <- TRUE
} else {
    thr <- suppressWarnings(as.numeric(hypermutator_max))
    if (is.na(thr) || thr <= 0) stop("hypermutator_max deve essere 'none' o un numero positivo, trovato: ", hypermutator_max)
    per_sample$kept <- per_sample$n_coding_events <= thr
}
message(sprintf("campioni: %d, esclusi come ipermutati: %d (soglia sugli eventi codificanti: %s)",
                nrow(per_sample), sum(!per_sample$kept), hypermutator_max))
events <- events[events$sample %in% per_sample$sample[per_sample$kept], ]

if (event_annotation == "reference") {
    # classe = classe del sito (posizione + base mutata) nei trascritti della modalita', stesse regole dei siti
    source("workflow/scripts/helpers/coding_sites_lib.R")
    g <- cs_read_transcripts(input_path_gtf, transcript_mode)
    G <- cs_read_genome(input_path_genome)
    events$sev <- 0L; n_ref_mismatch <- 0L
    for (ch in intersect(unique(cs_strip_chr(g$chr)), unique(events$chr[events$snv]))) {
        if (!ch %in% names(G$names)) next
        chseq <- G$seq[[G$names[[ch]]]]
        st <- cs_sites_chrom(g[cs_strip_chr(g$chr) == ch, ], chseq)
        ii <- which(events$chr == ch & events$snv)
        ref <- strsplit(as.character(chseq[events$pos[ii]]), "")[[1]]
        ok <- ref == events$wt[ii]; n_ref_mismatch <- n_ref_mismatch + sum(!ok)
        m <- match(paste(events$pos[ii], events$mt[ii]), paste(st$pos, st$alt))
        events$sev[ii] <- ifelse(ok & !is.na(m), st$sev[m] + 1L, 0L)        # sev della libreria 1..4 -> 2..5 (syn..trunc)
    }
    message(sprintf("riannotazione sulla sequenza di riferimento: SNV con base di riferimento COSMIC diversa dal genoma (esclusi): %d", n_ref_mismatch))
}
# classe per le colonne nuove: solo sostituzioni singole, e solo eventi con un trascritto della modalita'
events$cls <- ifelse(events$snv & events$sev >= 2L, CLASS_NAMES[pmax(events$sev, 1L)], NA_character_)
msg <- table(factor(ifelse(!events$snv, "non SNV (indel o multi-base)",
                    ifelse(events$sev == 0L, paste0("SNV non codificante nei trascritti ", transcript_mode),
                           paste0("SNV ", CLASS_NAMES[pmax(events$sev, 1L)])))))
message("eventi per categoria (modalita' ", transcript_mode, ", annotazione ", event_annotation, "):\n", paste(sprintf("  %-40s %d", names(msg), as.integer(msg)), collapse = "\n"))

# ---------------- assegnazione ai bin ----------------
backbone <- read.delim(file = input_path_backbone, header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)
backbone <- data.frame(chr = strip_chr(backbone$V1), start = backbone$V2 + 1, end = backbone$V3, binID = backbone$V4)

assigned <- list()
for (ch in intersect(unique(backbone$chr), unique(events$chr))) {
    b <- backbone[backbone$chr == ch, ]
    e <- events[events$chr == ch, ]
    hits <- IRanges::findOverlaps(IRanges::IRanges(start = e$pos, width = 1),
                                  IRanges::IRanges(start = b$start, end = b$end))
    if (length(hits) == 0) next
    assigned[[ch]] <- data.frame(binID = b$binID[S4Vectors::subjectHits(hits)],
                                 sample = e$sample[S4Vectors::queryHits(hits)],
                                 mut_id = e$mut_id[S4Vectors::queryHits(hits)],
                                 cls = e$cls[S4Vectors::queryHits(hits)],
                                 stringsAsFactors = FALSE)
}
assigned <- do.call(rbind, assigned)
if (is.null(assigned)) assigned <- data.frame(binID = character(), sample = character(), mut_id = character(), cls = character())

bf <- factor(assigned$binID, levels = unique(assigned$binID))
output_df <- data.frame(binID = levels(bf),
                        mut_n_events  = as.integer(table(bf)),
                        mut_n_samples = as.integer(tapply(assigned$sample, bf, function(x) length(unique(x)))),
                        mut_n_sites   = as.integer(tapply(assigned$mut_id, bf, function(x) length(unique(x)))),
                        stringsAsFactors = FALSE)
for (cl in c("syn", "mis", "trunc", "other_coding"))
    output_df[[paste0("mut_n_", cl)]] <- as.integer(tapply(assigned$cls %in% cl, bf, sum))
output_df <- output_df[order(output_df$binID, method = "radix"), ]                # ordine dei binID come dplyr::group_by (locale C)

write.table(x = output_df, file = output_path, quote = F, sep = "\t", row.names = F, col.names = T)
write.table(x = per_sample, file = output_path_samples, quote = F, sep = "\t", row.names = F, col.names = T)