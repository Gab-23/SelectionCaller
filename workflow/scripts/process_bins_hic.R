options(scipen=999)                                                               # unable scientific notation

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_backbone_raw                 <- args$input_path_backbone_raw
input_path_backbone_filtered            <- args$input_path_backbone_filtered
input_path_bins                         <- args$input_path_bins
output_path                             <- args$output_path

`%>%` <- magrittr::`%>%`


backbone_raw <- read.table(file = input_path_backbone_raw, header = F, sep = "\t")
backbone_filtered <- read.table(file = input_path_backbone_filtered, header = F, sep = "\t")
bins <- read.table(file = input_path_bins, header = F, sep = "\t")

colnames(backbone_raw) <- c("chr", "start", "end", "binID")
colnames(backbone_filtered) <- c("chr", "start", "end", "binID")
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
#
# NOTE (revised Option B design): matching is now done against backbone_raw (the pre-filter,
# continuous backbone that cload was actually run against - see compute_matrices_hic.smk) -
# this ALWAYS finds a match for every cooler bin, since backbone_raw has no gaps. A SEPARATE
# check against backbone_filtered's binID set then determines which of those bins are
# blacklisted (present in backbone_raw but absent from backbone_filtered, since
# filter_backbone.smk removes exactly the blacklisted regions) - these are set to binID = NA,
# same as before, but now derived by set-membership against the filtered backbone rather than
# a failed position match against it directly.
match_idx <- match(
    paste(bins$chr, bins$start, bins$end),
    paste(backbone_raw$chr, backbone_raw$start, backbone_raw$end)
)
stopifnot("some cooler bins could not be matched to the raw (pre-filter) backbone - this should never happen, since the raw backbone is continuous and cload was run against it" = all(!is.na(match_idx)))

raw_binID <- backbone_raw$binID[match_idx]
is_blacklisted <- !(raw_binID %in% backbone_filtered$binID)

bins_with_binID <- data.frame(
    chr = bins$chr,
    start = bins$start,
    end = bins$end,
    binID = ifelse(is_blacklisted, NA, raw_binID)
)

write.table(x = bins_with_binID, file = output_path, quote = F, col.names = F, row.names = F, sep = "\t", na = "NA")