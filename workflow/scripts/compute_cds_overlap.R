options(scipen=999)                                                               # unable scientific notation

# bp di sequenza codificante (CDS di geni protein_coding, GENCODE) per bin del backbone.
# Le CDS di trascritti diversi si sovrappongono: vengono fuse (IRanges::reduce) prima di
# contare, cosi' ogni bp codificante e' contato una volta sola.
# Coordinate: GTF 1-based chiuso; backbone BED 0-based semiaperto -> start + 1.
# I nomi dei cromosomi vengono confrontati senza il prefisso "chr" (GTF "chr1", backbone
# "chr1" o "1": funziona in entrambi i casi).

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_gtf          <- args$input_path_gtf
input_path_backbone     <- args$input_path_backbone
output_path             <- args$output_path

strip_chr <- function(x) sub("^chr", "", x)

gtf <- read.delim(file = gzfile(input_path_gtf), header = FALSE, sep = "\t", quote = "",
                  comment.char = "#", stringsAsFactors = FALSE,
                  colClasses = c("character", "NULL", "character", "integer", "integer",
                                 "NULL", "NULL", "NULL", "character"))
colnames(gtf) <- c("chr", "type", "start", "end", "attributes")
cds <- gtf[gtf$type == "CDS" & grepl('gene_type "protein_coding"', gtf$attributes, fixed = TRUE), c("chr", "start", "end")]
rm(gtf)
cds$chr <- strip_chr(cds$chr)

backbone <- read.delim(file = input_path_backbone, header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE)
backbone <- data.frame(chr = strip_chr(backbone$V1), start = backbone$V2 + 1, end = backbone$V3, binID = backbone$V4)

cds_bp <- setNames(numeric(nrow(backbone)), backbone$binID)
for (ch in intersect(unique(backbone$chr), unique(cds$chr))) {
    b <- backbone[backbone$chr == ch, ]
    bins <- IRanges::IRanges(start = b$start, end = b$end)
    exons <- IRanges::reduce(IRanges::IRanges(start = cds$start[cds$chr == ch], end = cds$end[cds$chr == ch]))
    hits <- IRanges::findOverlaps(bins, exons)
    if (length(hits) == 0) next
    w <- IRanges::width(IRanges::pintersect(bins[S4Vectors::queryHits(hits)], exons[S4Vectors::subjectHits(hits)]))
    per_bin <- tapply(w, S4Vectors::queryHits(hits), sum)
    cds_bp[b$binID[as.integer(names(per_bin))]] <- as.numeric(per_bin)
}

write.table(x = data.frame(binID = names(cds_bp), cds_bp = as.numeric(cds_bp)),
            file = output_path, quote = F, sep = "\t", row.names = F, col.names = T)