options(scipen=999)                                                               # unable scientific notation

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_backbone                     <- args$input_path_backbone
input_path_bins                         <- args$input_path_bins
output_path                             <- args$output_path

`%>%` <- magrittr::`%>%`


backbone <- read.table(file = input_path_backbone, header = F, sep = "\t")
bins <- read.table(file = input_path_bins, header = F, sep = "\t")

colnames(backbone) <- c("chr", "start", "end", "binID")
colnames(bins) <- c("chr", "start", "end", "unk")

# NOTE (fix): the original used merge(x=backbone, y=bins, by=c("chr","start","end")), which
# by default SORTS the result alphabetically by the "by" columns - "chr10" sorts before "chr2"
# alphabetically, even though it comes after it in every real assembly's actual chromosome
# order. Confirmed with a toy test using chr1/chr2/chr10: merge() reordered the rows so chr10
# preceded chr2. Since `bins` is cooler's own bin table - whose ROW ORDER *is* the contact
# matrix's row/column order - and process_hic.py assumes bins[i] corresponds to matrix row i,
# this reordering would silently misassign every downstream eigenvector to the wrong bin for
# any real genome with 10+ numbered chromosomes (i.e. every real human/mouse assembly).
# Fixed with an explicit, order-preserving lookup that keeps `bins`' original row order.
match_idx <- match(
    paste(bins$chr, bins$start, bins$end),
    paste(backbone$chr, backbone$start, backbone$end)
)
stopifnot("some cooler bins could not be matched to the backbone" = all(!is.na(match_idx)))

bins_with_binID <- data.frame(
    chr = bins$chr,
    start = bins$start,
    end = bins$end,
    binID = backbone$binID[match_idx]
)

write.table(x = bins_with_binID, file = output_path, quote = F, col.names = F, row.names = F, sep = "\t")
