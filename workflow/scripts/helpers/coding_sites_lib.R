# Libreria condivisa (27-28/10/2026): enumerazione delle SOSTITUZIONI SINGOLE POSSIBILI nei trascritti codificanti e loro
# classe, dalla sequenza di riferimento. Usata da:
#   compute_coding_sites.R             -> conta i siti per bin (esposizione delle classi)
#   compute_coding_mutations_counts.R  -> classifica gli EVENTI osservati con le stesse regole (riannotazione)
# cosi' eventi ed esposizione sono definiti nello stesso modo per costruzione.
# Regole (invariate rispetto a compute_coding_sites.R del 24/10/2026):
#   ORF = CDS + stop_codon del GTF; codone interno: sinonima / missense / nonsense (-> troncante);
#   primo codone e codone di stop: "altra codificante" (start_lost, stop_lost, stop_retained);
#   siti di splicing essenziali: 2 basi introniche a ogni confine tra esoni consecutivi (introni >= 4 bp) -> troncanti;
#   trascritti esclusi: CDS incompleta (tag cds_start_NF / cds_end_NF), selenoproteine (tag seleno), ORF non multipla
#   di 3 o < 6 bp, basi non ACGT, cromosoma assente dal genoma;
#   stessa sostituzione (posizione + base alternativa) in piu' trascritti / geni -> classe PIU' GRAVE, contata una volta.
# Gravita' (sev): 1 sinonima, 2 altra codificante, 3 missense, 4 troncante.
suppressMessages({ library(Biostrings) })
CS_CLASS_NAMES <- c("syn", "other_coding", "mis", "trunc")
CS_BASES <- c("A", "C", "G", "T")
cs_strip_chr <- function(x) sub("^chr", "", x)
.cs_comp <- c(A = "T", C = "G", G = "C", T = "A")

.cs_lut <- local({
    codons <- names(GENETIC_CODE); aa <- GENETIC_CODE
    classify <- function(ref, alt, context) {
        if (context == "first") return(2L)
        if (context == "last")  return(2L)
        a1 <- aa[[ref]]; a2 <- aa[[alt]]
        if (a1 == "*") return(2L)
        if (a2 == "*") return(4L)
        if (a1 == a2) 1L else 3L
    }
    LUT <- array(NA_integer_, dim = c(64, 3, 4, 3), dimnames = list(codons, 1:3, CS_BASES, c("internal", "first", "last")))
    for (ctx in c("internal", "first", "last")) for (cd in codons) for (p in 1:3) for (b in CS_BASES) {
        ch <- strsplit(cd, "")[[1]]
        if (ch[p] == b) next
        alt <- ch; alt[p] <- b
        LUT[cd, p, b, ctx] <- classify(cd, paste(alt, collapse = ""), ctx)
    }
    list(LUT = LUT, codons = codons)
})

cs_read_transcripts <- function(gtf_path, transcript_mode) {
    # Righe del GTF (transcript / exon / CDS / stop_codon) dei trascritti usati, dopo le esclusioni per tag.
    if (!transcript_mode %in% c("canonical", "all")) stop("transcript_mode deve essere 'canonical' o 'all', trovato: ", transcript_mode)
    g <- read.delim(gzfile(gtf_path), header = FALSE, comment.char = "#", quote = "", sep = "\t",
                    colClasses = c("character", "NULL", "character", "integer", "integer", "NULL", "character", "NULL", "character"),
                    col.names = c("chr", "src", "type", "start", "end", "score", "strand", "frame", "attr"))
    g <- g[g$type %in% c("transcript", "exon", "CDS", "stop_codon"), ]
    g$tx <- sub('.*transcript_id "([^".]+)[^"]*".*', "\\1", g$attr)
    tr <- g[g$type == "transcript", ]
    has_tag <- function(a, t) grepl(paste0('tag "', t, '"'), a, fixed = TRUE)
    tr$canonical <- has_tag(tr$attr, "Ensembl_canonical")
    tr$incomplete <- has_tag(tr$attr, "cds_start_NF") | has_tag(tr$attr, "cds_end_NF")
    tr$seleno <- has_tag(tr$attr, "seleno")
    cds_tx <- unique(g$tx[g$type == "CDS"])
    tr <- tr[tr$tx %in% cds_tx, ]
    if (transcript_mode == "canonical") tr <- tr[tr$canonical, ]
    n0 <- nrow(tr)
    tr <- tr[!tr$incomplete & !tr$seleno, ]
    message(sprintf("trascritti con CDS nella modalita' %s: %d; esclusi (CDS incompleta per tag o selenoproteina): %d", transcript_mode, n0, n0 - nrow(tr)))
    g[g$tx %in% tr$tx, c("chr", "type", "start", "end", "strand", "tx")]
}

cs_read_genome <- function(genome_path) {
    genome <- readDNAStringSet(genome_path)
    names(genome) <- sub("\\s.*$", "", names(genome))
    list(seq = genome, names = setNames(names(genome), cs_strip_chr(names(genome))))
}

cs_sites_chrom <- function(gc, chseq) {
    # Sostituzioni possibili dei trascritti di UN cromosoma: data.frame(pos, alt = base sul filamento +, sev), una riga per
    # (pos, alt) con la classe piu' grave tra i trascritti; attributo "n_bad" = trascritti esclusi (ORF o basi).
    LUT <- .cs_lut$LUT; codons <- .cs_lut$codons; comp <- .cs_comp
    by_tx <- split(gc, gc$tx)
    out_pos <- vector("list", length(by_tx)); out_alt <- out_pos; out_sev <- out_pos; n_bad <- 0L
    for (k in seq_along(by_tx)) {
        d <- by_tx[[k]]; strand <- d$strand[1]
        orf <- d[d$type %in% c("CDS", "stop_codon"), ]; orf <- orf[order(orf$start), ]
        pos <- unlist(mapply(seq, orf$start, orf$end, SIMPLIFY = FALSE), use.names = FALSE)
        pos <- unique(pos)
        if (length(pos) %% 3 != 0 || length(pos) < 6) { n_bad <- n_bad + 1L; next }
        if (strand == "-") pos <- rev(pos)
        s <- strsplit(as.character(chseq[sort(pos)]), "")[[1]]
        s <- s[match(pos, sort(pos))]
        if (strand == "-") s <- unname(comp[s])
        if (any(is.na(s))) { n_bad <- n_bad + 1L; next }
        L <- length(s); ncod <- L / 3
        cod <- paste0(s[seq(1, L, 3)], s[seq(2, L, 3)], s[seq(3, L, 3)])
        ctx <- rep("internal", ncod); ctx[1] <- "first"; ctx[ncod] <- "last"
        ci <- (seq_len(L) - 1) %/% 3 + 1; pin <- (seq_len(L) - 1) %% 3 + 1
        idx <- rep(seq_len(L), each = 4)
        Bg <- rep(CS_BASES, L)
        Bt <- if (strand == "-") unname(comp[Bg]) else Bg
        keep <- Bt != s[idx]
        sev <- LUT[cbind(match(cod[ci[idx]], codons), pin[idx], match(Bt, CS_BASES), match(ctx[ci[idx]], c("internal", "first", "last")))]
        ex <- d[d$type == "exon", ]; ex <- ex[order(ex$start), ]
        sp <- integer(0)
        if (nrow(ex) > 1) for (j in 1:(nrow(ex) - 1)) {
            i1 <- ex$end[j] + 1; i2 <- ex$start[j + 1] - 1
            if (i2 - i1 + 1 >= 4) sp <- c(sp, i1, i1 + 1, i2 - 1, i2)
        }
        refsp <- if (length(sp)) strsplit(as.character(chseq[sp]), "")[[1]] else character(0)
        Psp <- rep(sp, each = 4); Bsp <- rep(CS_BASES, length(sp)); ksp <- Bsp != rep(refsp, each = 4)
        out_pos[[k]] <- c(pos[idx][keep], Psp[ksp])
        out_alt[[k]] <- c(Bg[keep], Bsp[ksp])
        out_sev[[k]] <- c(sev[keep], rep(4L, sum(ksp)))
    }
    P <- unlist(out_pos); A <- unlist(out_alt); S <- unlist(out_sev)
    if (!length(P)) { r <- data.frame(pos = numeric(0), alt = character(0), sev = integer(0)); attr(r, "n_bad") <- n_bad; return(r) }
    key <- P * 4 + match(A, CS_BASES)
    o <- order(key, -S); key <- key[o]; S <- S[o]; P <- P[o]; A <- A[o]
    first <- !duplicated(key)
    r <- data.frame(pos = P[first], alt = A[first], sev = S[first], stringsAsFactors = FALSE); attr(r, "n_bad") <- n_bad
    r
}