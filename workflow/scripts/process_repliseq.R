source("workflow/scripts/helpers/CPM.R")
source("workflow/scripts/helpers/find_contiguous_subsets.R")

options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_dir                           <- args$input_path_dir
input_path_samples                       <- args$input_path_samples
input_path_library_size                  <- args$input_path_library_size
input_path_backbone_mask                 <- args$input_path_backbone_mask
output_path_data                         <- args$output_path_data
output_path_plot_distribution            <- args$output_path_plot_distribution
output_path_plot_landscape               <- args$output_path_plot_landscape
output_path_plot_raw_smooth              <- args$output_path_plot_raw_smooth
window_size                              <- as.numeric(args$window_size)
smoothing_size                           <- as.numeric(args$smoothing_size)
assembly_code                            <- args$assembly_code
window_mode                              <- args$window_mode
step                                     <- args$step
mapq_thr                                 <- args$mapq_thr
cov_thr                                  <- args$cov_thr
chr_to_plot                              <- 17

`%>%` <- magrittr::`%>%`

files <- list.files(input_path_dir, recursive = T, full.names = T, pattern = "filtered")
# NOTE (fix): the coverage files under input_path_dir accumulate across EVERY
# assembly_code/window_mode/window_size/step/mapq_thr/cov_thr combination ever run (e.g. both
# window_size=100000 and window_size=20000, if both have been run) - filtering only by the
# literal word "filtered" (present in every such filename regardless of parameters) let files
# from a DIFFERENT parameter combination slip in, causing "expected exactly 1 coverage file...
# found 2" once more than one combination coexisted on disk. Confirmed directly on a real run.
# Narrowed to the EXACT filename this specific run expects (anchored at the end of the path,
# not just a loose substring - avoids e.g. window_size=20 matching inside window_size=200).
expected_suffix <- paste0(assembly_code, "_", window_mode, "_", window_size, "_", step, "_", mapq_thr, "_", cov_thr, ".filtered.bed")
files <- files[grepl(paste0(expected_suffix, "$"), files)]
samples <- read.table(file = input_path_samples, header = T, sep = "\t")
library_size_table <- read.table(file = input_path_library_size, header = T, sep = "\t")

mask <- read.table(file = input_path_backbone_mask, header = F, sep = "\t")
colnames(mask) <- c("chr", "start", "end", "binID")

samples <- merge(samples, library_size_table, by = "SRA_id")

sample_types <- c("ES", "LS", "G1")

binIDs_list <- list()
CPM_list <- list()

for (idx in seq(length(sample_types))){
    
    sample_type <- sample_types[[idx]]

    pattern <- samples[samples$sample_type == sample_type, "SRA_id"]
    library_size <- samples[samples$sample_type == sample_type, "library_size"]

    file <- grep(pattern = paste(pattern, collapse = "|"), x = files, value = T)   # fix: pattern can have >1 SRA_id (multiple replicates per sample_type) - grep() previously silently used only the first element of a vector pattern, dropping all other replicate files with just a warning, no error
    if (length(file) != 1){                                                       # fail loudly instead of guessing how to merge multiple replicate coverage
        stop(sprintf(                                                            # files (sum? average CPM?) - that is a scientific design decision this
            "process_repliseq.R: expected exactly 1 coverage file for sample_type '%s', found %d: %s",
            sample_type, length(file), paste(file, collapse = ", ")
        ))
    }
    
    coverage <- read.table(file, header = F, sep = "\t")
    colnames(coverage) <- c("chr", "start", "end", "binID", "coverage")

    norm_coverage <- CPM(coverage, library_size)
    
    binIDs_list[[sample_type]] <- norm_coverage$binID
    CPM_list[[sample_type]] <- norm_coverage

}

common_binIDs <- Reduce(f = intersect, x = binIDs_list)

ES <- CPM_list[["ES"]]
LS <- CPM_list[["LS"]]
G1 <- CPM_list[["G1"]]

ES.f <- ES[ES$binID %in% common_binIDs, ]
LS.f <- LS[LS$binID %in% common_binIDs, ]
G1.f <- G1[G1$binID %in% common_binIDs, ]

stopifnot(all(ES.f$binID == G1.f$binID))
stopifnot(all(LS.f$binID == G1.f$binID))

stopifnot(all(ES.f$chr == G1.f$chr) & all(ES.f$start == G1.f$start) & all(ES.f$end == G1.f$end))
stopifnot(all(LS.f$chr == G1.f$chr) & all(LS.f$start == G1.f$start) & all(LS.f$end == G1.f$end))

ES_G1_norm_coverage <- ES.f$CPM / G1.f$CPM
LS_G1_norm_coverage <- LS.f$CPM / G1.f$CPM

EL_coverage <- ES_G1_norm_coverage / LS_G1_norm_coverage

log2_EL_coverage <- log2(EL_coverage)

replication_timing <- data.frame(binID = common_binIDs,
                                 start = G1.f$start,
                                 end   = G1.f$end,
                                 replication.timing = log2_EL_coverage)

replication_timing <- tidyr::separate_wider_delim(data = replication_timing,
                                                  cols = binID,
                                                  delim = "_", 
                                                  names = c("chr", "pos"))

replication_timing$chr <- paste0("chr", replication_timing$chr)                    # fix: binID's chr-part has the "chr" prefix stripped (see
                                                                                    # generate_backbone.R: binID <- paste0(gsub("chr","",chr), ...)),
                                                                                    # but the backbone MASK file (read a few lines above) uses real
                                                                                    # chromosome names ("chr1", not "1"). Without re-adding the
                                                                                    # prefix here, the GRanges built from replication_timing below
                                                                                    # would have seqnames "1".."22" while mask_gr has "chr1".."chr22" -
                                                                                    # findOverlaps() would then ALWAYS return zero hits (confirmed with
                                                                                    # a toy test: the final output was silently empty, header only).

replication_timing$replication.timing.smoothed <- NA

for (chr in unique(replication_timing$chr)){
  
  chr_mask <- replication_timing$chr == chr
  chr_subset <- replication_timing[chr_mask, , drop = FALSE]

  bins_in_smoothing_size <- (smoothing_size / window_size)                        # how many bins of size window_size to obtain smoothing_size

  contiguous_subsets <- find_contiguous_subsets(chr_subset)

  for (contiguous_subset in contiguous_subsets){
    contiguous_subset_mask <- ((replication_timing$chr == chr) & (replication_timing$pos %in% contiguous_subset$pos))

    bins_per_contiguous_subset <- nrow(contiguous_subset)

    span = min(1,(bins_in_smoothing_size / bins_per_contiguous_subset))           # if bins_per_contiguous subset < bins_in_smoothing_size, i.e. there is not enough bins to cover smoothing size, smooth for all the dataset 

    if (bins_per_contiguous_subset > 5){                                          # if there is "enough" bins for smoothing, smooth, else avoid

    smoothener <- loess(replication.timing ~ pos,
                      data = contiguous_subset,
                      span = span,
                      degree = 1)

    replication_timing$replication.timing.smoothed[contiguous_subset_mask] <- predict(smoothener)
    } else {
    replication_timing$replication.timing.smoothed[contiguous_subset_mask] <- replication_timing$replication.timing[contiguous_subset_mask]
    }
  }
}

replication_timing$binID <- paste0(replication_timing$chr, "_", replication_timing$pos)

hist <- ggplot2::ggplot(replication_timing, ggplot2::aes(x = replication.timing)) +
  ggplot2::geom_histogram(
    bins = 50,
    fill = "skyblue",
    color = "black",
    alpha = 0.7
  ) +
  ggplot2::labs(
    title = "log2(E/L) ratio",
    x = "log2(E/L)",
    y = "Frequency"
  ) +
  ggplot2::scale_x_continuous(n.breaks = 30) +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
  )  

landscape <- ggplot2::ggplot(
  replication_timing %>%
    dplyr::mutate(
      chr = sub("_.*", "", binID),
      chr_num = as.numeric(sub("chr", "", chr)),
      pos = as.numeric(sub(".*_", "", binID))
    ) %>%
    dplyr::arrange(chr_num, pos) %>%
    dplyr::mutate(chr = factor(chr, levels = unique(chr))),
  ggplot2::aes(x = pos, y = replication.timing.smoothed)
  ) +
  ggplot2::geom_line() +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.2
  ) +
  ggplot2::facet_wrap(~ chr, scales = "free_x") +
  ggplot2::labs(
    x = "Genomic position",
    y = "Smoothed Replication Timing (log2)",
    title = "Replication timing across chromosomes"
  ) +
  ggplot2::theme_bw()

raw_smooth <- ggplot2::ggplot(
  replication_timing %>%
    dplyr::mutate(
      chr = sub("_.*", "", binID),
      chr_num = as.numeric(sub("chr", "", chr)),
      pos = as.numeric(sub(".*_", "", binID))
    ) %>%
    dplyr::filter(chr_num == chr_to_plot) %>%
    dplyr::arrange(chr_num, pos) %>%
    dplyr::mutate(chr = factor(chr, levels = unique(chr))),
  ) +
  ggplot2::geom_point(ggplot2::aes(x = pos, y = replication.timing), color = "lightblue", size = 1) +
  ggplot2::geom_line(ggplot2::aes(x = pos, y = replication.timing.smoothed), color = "red", linewidth = 1) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.2
  ) +
  ggplot2::theme_classic() +
  ggplot2::labs(y = "log2(Early/Late)", 
                title = paste0("Raw vs Smoothed signal [chr",chr_to_plot, "]"))

ggplot2::ggsave(filename = output_path_plot_distribution, plot = hist, width = 8, height = 6, dpi = 300)
ggplot2::ggsave(filename = output_path_plot_landscape, plot = landscape, width = 24, height = 10, units = "in")
ggplot2::ggsave(filename = output_path_plot_raw_smooth, plot = raw_smooth, width = 24, height = 10, units = "in")

replication_timing_gr <- GenomicRanges::makeGRangesFromDataFrame(
  replication_timing,
  seqnames.field = "chr",
  start.field = "start",
  end.field = "end",
  starts.in.df.are.0based = TRUE,
  keep.extra.columns = TRUE
)

mask_gr <- GenomicRanges::makeGRangesFromDataFrame(
  mask,
  seqnames.field = "chr",
  start.field = "start",
  end.field = "end",
  starts.in.df.are.0based = TRUE,
  keep.extra.columns = TRUE
)

hits <- GenomicRanges::findOverlaps(query = mask_gr, subject = replication_timing_gr, type = "any")

overlap_table <- data.frame(
  mask_idx     = S4Vectors::queryHits(hits),
  backbone_idx = S4Vectors::subjectHits(hits)
)
overlap_table$replication.timing.smoothed <- replication_timing_gr$replication.timing.smoothed[overlap_table$backbone_idx]

averaged <- overlap_table %>%
  dplyr::select(-backbone_idx) %>%
  dplyr::group_by(mask_idx) %>%
  dplyr::summarise(avg_RT = mean(replication.timing.smoothed))

output_data <- data.frame(
  binID                       = mask$binID[averaged$mask_idx],
  replication.timing.smoothed = averaged$avg_RT
)

write.table(x = output_data,
            file = output_path_data, 
            quote = F, 
            sep = "\t", 
            row.names = F, 
            col.names = T)