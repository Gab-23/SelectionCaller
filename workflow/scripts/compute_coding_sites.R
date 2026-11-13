options(scipen=999)                                                               # unable scientific notation

# Esposizione per classe di mutazione: per ogni bin, il numero di SOSTITUZIONI SINGOLE POSSIBILI che
# produrrebbero ciascuna classe, contate dalla sequenza di riferimento. E' il denominatore (offset) dei
# conteggi per classe di compute_coding_mutations_counts.R, con le stesse definizioni:
#   sites_syn            sinonime (stesso amminoacido)
#   sites_mis            missense (amminoacido diverso)
#   sites_trunc          troncanti: nonsense (nasce uno stop) + siti di splicing essenziali (le 2 basi
#                        introniche adiacenti a ogni confine esone-introne: donatore +1/+2, accettore -1/-2,
#                        come splice_donor_variant / splice_acceptor_variant di COSMIC/VEP)
#   sites_other_coding   primo codone (start_lost), codone di stop (stop_lost / stop_retained_variant)
# Trascritti (transcript_mode): canonical = tag Ensembl_canonical; all = tutti i trascritti con CDS (esclusioni nel log).
# Sovrapposizioni: una sostituzione possibile presente in piu' trascritti o geni prende la classe PIU' GRAVE e conta UNA volta.
# NOTE (28/10/2026): l'enumerazione e' stata spostata in helpers/coding_sites_lib.R, condivisa con
# compute_coding_mutations_counts.R (che ora classifica gli eventi con le stesse regole). Il risultato e' identico.
# Output: una riga per OGNI bin del backbone (0 dove non ci sono siti).

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args
input_path_genome   <- args$input_path_genome     # FASTA del genoma di riferimento
input_path_gtf      <- args$input_path_gtf        # GENCODE GTF (gz)
input_path_backbone <- args$input_path_backbone
transcript_mode     <- if (is.null(args$transcript_mode)) "canonical" else as.character(args$transcript_mode)
output_path         <- args$output_path
source("workflow/scripts/helpers/coding_sites_lib.R")

g <- cs_read_transcripts(input_path_gtf, transcript_mode)
G <- cs_read_genome(input_path_genome)

sites <- list(); n_bad_len <- 0L; n_no_seq <- 0L
for (ch in unique(cs_strip_chr(g$chr))) {
    gc <- g[cs_strip_chr(g$chr) == ch, ]
    if (!ch %in% names(G$names)) { n_no_seq <- n_no_seq + length(unique(gc$tx)); next }
    s <- cs_sites_chrom(gc, G$seq[[G$names[[ch]]]])
    n_bad_len <- n_bad_len + attr(s, "n_bad")
    if (nrow(s)) sites[[ch]] <- data.frame(chr = ch, pos = s$pos, sev = s$sev)
}
message(sprintf("trascritti esclusi: lunghezza ORF non multipla di 3 o basi non ACGT = %d; cromosoma assente dal genoma = %d", n_bad_len, n_no_seq))
sites <- do.call(rbind, sites)

backbone <- read.delim(file = input_path_backbone, header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)
backbone <- data.frame(chr = cs_strip_chr(backbone$V1), start = backbone$V2 + 1, end = backbone$V3, binID = backbone$V4)
res <- data.frame(binID = backbone$binID, sites_syn = 0L, sites_mis = 0L, sites_trunc = 0L, sites_other_coding = 0L)
for (ch in intersect(unique(backbone$chr), unique(sites$chr))) {
    bi <- which(backbone$chr == ch); b <- backbone[bi, ]; s <- sites[sites$chr == ch, ]
    hits <- IRanges::findOverlaps(IRanges::IRanges(start = s$pos, width = 1), IRanges::IRanges(start = b$start, end = b$end))
    if (length(hits) == 0) next
    tab <- table(factor(bi[S4Vectors::subjectHits(hits)], levels = bi), factor(CS_CLASS_NAMES[s$sev[S4Vectors::queryHits(hits)]], levels = CS_CLASS_NAMES))
    res[bi, paste0("sites_", CS_CLASS_NAMES)] <- res[bi, paste0("sites_", CS_CLASS_NAMES)] + matrix(as.integer(tab), ncol = 4)
}
message(sprintf("siti totali: syn %d, mis %d, trunc %d, other_coding %d", sum(res$sites_syn), sum(res$sites_mis), sum(res$sites_trunc), sum(res$sites_other_coding)))
write.table(x = res, file = output_path, quote = F, sep = "\t", row.names = F, col.names = T)