source("workflow/scripts/helpers/CPM.R")
source("workflow/scripts/helpers/find_contiguous_subsets.R")

options(scipen=999)                                                               # unable scientific notation           

# NOTE (30/10/2026): REPLICATI e LINEA CELLULARE.
# - Il file dei campioni (uno per linea cellulare) ha le colonne SRA_id, sample_type, cell_line e, facoltativa, replicate.
# - COMBINAZIONE dei replicati, per sample_type (ES, LS, G1): si SOMMANO le coperture dei replicati, bin per bin, e le loro
#   dimensioni di libreria, PRIMA della CPM (equivale a mettere insieme le letture). Si usano i bin presenti in TUTTI i
#   replicati (dopo il filtro di copertura di filter_coverage_repliseq), come gia' si faceva tra sample_type.
#   Con un solo campione per sample_type il risultato e' IDENTICO alla versione precedente.
# - Il file di uscita e i tre grafici di sempre usano la combinazione. In piu', se ci sono >= 2 replicati, il calcolo viene
#   ripetuto per ogni replicato e un PDF aggiuntivo (output_path_plot_replicates) mostra, in griglia: andamento lisciato
#   per cromosoma (righe = replicati e combinazione), distribuzioni, confronto tra replicati (correlazione) e segnale grezzo
#   contro lisciato sul cromosoma chr_to_plot. Con un solo replicato il PDF contiene una pagina che lo dice.
# Il calcolo di log2(E/L), la lisciatura e la media sui bin della maschera sono INVARIATI.

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_dir                           <- args$input_path_dir
input_path_samples                       <- args$input_path_samples
input_path_library_size                  <- args$input_path_library_size
input_path_backbone_mask                 <- args$input_path_backbone_mask
output_path_data                         <- args$output_path_data
output_path_plot_distribution            <- args$output_path_plot_distribution
output_path_plot_landscape               <- args$output_path_plot_landscape
output_path_plot_raw_smooth              <- args$output_path_plot_raw_smooth
output_path_plot_replicates              <- args$output_path_plot_replicates
window_size                              <- as.numeric(args$window_size)
smoothing_size                           <- as.numeric(args$smoothing_size)
assembly_code                            <- args$assembly_code
window_mode                              <- args$window_mode
step                                     <- args$step
mapq_thr                                 <- args$mapq_thr
cov_thr                                  <- args$cov_thr
chr_to_plot                              <- 17

`%>%` <- magrittr::`%>%`

bins_in_smoothing_size <- (smoothing_size / window_size)                        # how many bins of size window_size to obtain smoothing_size
  min_bins_in_smoothing_size <- 4
  if (bins_in_smoothing_size < min_bins_in_smoothing_size) {
    stop(sprintf(paste0(
      "smoothing_size (%s bp) covers only %.2f bin(s) of window_size (%s bp): loess (degree = 1) ",
      "needs at least %d bins per local fit (fewer either fails with 'span is too small' or ",
      "gives degenerate fits). Use smoothing_size >= %s bp for this window_size."),
      format(smoothing_size, scientific = FALSE), bins_in_smoothing_size,
      format(window_size, scientific = FALSE), min_bins_in_smoothing_size,
      format(min_bins_in_smoothing_size * window_size, scientific = FALSE)))
  }

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
if (!"replicate" %in% colnames(samples)) samples$replicate <- 1                  # colonna facoltativa: senza, un solo replicato
library_size_table <- read.table(file = input_path_library_size, header = T, sep = "\t")

mask <- read.table(file = input_path_backbone_mask, header = F, sep = "\t")
colnames(mask) <- c("chr", "start", "end", "binID")

n_before <- nrow(samples)
samples <- merge(samples, library_size_table, by = "SRA_id")
if (nrow(samples) != n_before) stop("process_repliseq.R: alcuni SRA_id del file dei campioni non hanno una dimensione di libreria")

sample_types <- c("ES", "LS", "G1")

coverage_of <- function(sra){
    # file di copertura di UN campione: esattamente uno
    f <- files[grepl(sra, files, fixed = TRUE)]
    if (length(f) != 1) stop(sprintf("process_repliseq.R: expected exactly 1 coverage file for SRA_id '%s', found %d: %s",
                                     sra, length(f), paste(f, collapse = ", ")))
    cov <- read.table(f, header = F, sep = "\t")
    colnames(cov) <- c("chr", "start", "end", "binID", "coverage")
    cov
}

pooled_coverage <- function(s){
    # coperture SOMMATE dei campioni s (stesso sample_type) sui bin comuni a tutti; dimensione di libreria sommata
    covs <- lapply(s$SRA_id, coverage_of)
    common <- Reduce(intersect, lapply(covs, function(x) x$binID))
    base <- covs[[1]][covs[[1]]$binID %in% common, c("chr", "start", "end", "binID")]
    base$coverage <- 0
    for (x in covs) base$coverage <- base$coverage + x$coverage[match(base$binID, x$binID)]
    list(coverage = base, library_size = sum(s$library_size))
}

compute_replication_timing <- function(s){
    # log2(E/L) e lisciatura, per l'insieme di campioni s (uno o piu' replicati): stesso calcolo di prima
    binIDs_list <- list()
    CPM_list <- list()
    for (sample_type in sample_types){
        st <- s[s$sample_type == sample_type, ]
        if (nrow(st) == 0) stop(sprintf("process_repliseq.R: nessun campione per sample_type '%s'", sample_type))
        pooled <- pooled_coverage(st)
        norm_coverage <- CPM(pooled$coverage, pooled$library_size)
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
                                                                                        # but the backbone MASK file uses real chromosome names ("chr1",
                                                                                        # not "1"); without the prefix findOverlaps() returns zero hits.

    replication_timing$replication.timing.smoothed <- NA

    for (chr in unique(replication_timing$chr)){
      
      chr_mask <- replication_timing$chr == chr
      chr_subset <- replication_timing[chr_mask, , drop = FALSE]

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
    replication_timing
}

# ---------------- combinazione di tutti i replicati (uscita principale) ----------------
replication_timing <- compute_replication_timing(samples)

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

order_chr <- function(df) {
  df %>%
    dplyr::mutate(
      chr = sub("_.*", "", binID),
      chr_num = suppressWarnings(as.numeric(sub("chr", "", chr))),
      chr_num = ifelse(is.na(chr_num), ifelse(chr == "chrX", 23, 24), chr_num),
      pos = as.numeric(sub(".*_", "", binID))
    ) %>%
    dplyr::arrange(chr_num, pos) %>%
    dplyr::mutate(chr = factor(chr, levels = unique(chr)))
}

landscape <- ggplot2::ggplot(
  order_chr(replication_timing),
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
  order_chr(replication_timing) %>% dplyr::filter(chr_num == chr_to_plot)
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

# ---------------- replicati: stesso calcolo per ciascuno, grafici in griglia ----------------
replicates <- sort(unique(samples$replicate))
pdf(output_path_plot_replicates, width = 24, height = 4 + 3 * (length(replicates) + 1))
if (length(replicates) < 2) {
  print(ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0, size = 8,
        label = paste0("Un solo replicato (", paste(unique(samples$cell_line), collapse = ", "), "): nessun confronto tra replicati")) +
        ggplot2::theme_void())
} else {
  per_rep <- lapply(replicates, function(r) compute_replication_timing(samples[samples$replicate == r, ]))
  sets <- c(paste0("replicate ", replicates), "combined")
  all_sets <- dplyr::bind_rows(c(per_rep, list(replication_timing)), .id = "set_idx") %>%
    dplyr::mutate(set = factor(sets[as.integer(set_idx)], levels = sets))
  all_sets <- order_chr(all_sets)
  # 1. andamento lisciato per cromosoma: righe = replicati e combinazione
  print(ggplot2::ggplot(all_sets, ggplot2::aes(x = pos, y = replication.timing.smoothed)) +
          ggplot2::geom_line(linewidth = 0.3) +
          ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.2) +
          ggplot2::facet_grid(set ~ chr, scales = "free_x", space = "free_x") +
          ggplot2::labs(x = "Genomic position", y = "Smoothed log2(E/L)",
                        title = "Replication timing per replicate and combined") +
          ggplot2::theme_bw() + ggplot2::theme(axis.text.x = ggplot2::element_blank()))
  # 2. distribuzioni
  print(ggplot2::ggplot(all_sets, ggplot2::aes(x = replication.timing)) +
          ggplot2::geom_histogram(bins = 50, fill = "skyblue", color = "black", alpha = 0.7) +
          ggplot2::facet_wrap(~ set, nrow = 1) +
          ggplot2::labs(x = "log2(E/L)", y = "Frequency", title = "log2(E/L) per replicate and combined") +
          ggplot2::theme_minimal(base_size = 14))
  # 3. confronto tra replicati (segnale lisciato, bin comuni)
  wide <- Reduce(function(a, b) merge(a, b, by = "binID"),
                 lapply(seq_along(per_rep), function(i) {
                   d <- per_rep[[i]][, c("binID", "replication.timing.smoothed")]; colnames(d)[2] <- paste0("rep", replicates[i]); d }))
  pairs <- utils::combn(paste0("rep", replicates), 2, simplify = FALSE)
  pair_df <- dplyr::bind_rows(lapply(pairs, function(p) data.frame(pair = paste(p, collapse = " vs "), x = wide[[p[1]]], y = wide[[p[2]]])))
  lab <- pair_df %>% dplyr::group_by(pair) %>%
    dplyr::summarise(label = sprintf("Pearson %.3f, Spearman %.3f, n = %d", cor(x, y), cor(x, y, method = "spearman"), dplyr::n()), .groups = "drop")
  message("concordanza tra replicati (segnale lisciato):\n", paste(lab$pair, lab$label, sep = ": ", collapse = "\n"))
  print(ggplot2::ggplot(pair_df, ggplot2::aes(x = x, y = y)) +
          ggplot2::geom_point(alpha = 0.2, size = 0.5) + ggplot2::geom_abline(linetype = "dashed", color = "red") +
          ggplot2::geom_text(data = lab, ggplot2::aes(x = -Inf, y = Inf, label = label), hjust = -0.05, vjust = 1.5, inherit.aes = FALSE) +
          ggplot2::facet_wrap(~ pair) + ggplot2::coord_equal() +
          ggplot2::labs(x = "smoothed log2(E/L)", y = "smoothed log2(E/L)", title = "Agreement between replicates") +
          ggplot2::theme_bw())
  # 4. grezzo contro lisciato su chr_to_plot
  rs <- all_sets %>% dplyr::filter(chr_num == chr_to_plot)
  if (nrow(rs) > 0) print(ggplot2::ggplot(rs) +
          ggplot2::geom_point(ggplot2::aes(x = pos, y = replication.timing), color = "lightblue", size = 1) +
          ggplot2::geom_line(ggplot2::aes(x = pos, y = replication.timing.smoothed), color = "red", linewidth = 1) +
          ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.2) +
          ggplot2::facet_grid(set ~ .) + ggplot2::theme_classic() +
          ggplot2::labs(y = "log2(Early/Late)", title = paste0("Raw vs Smoothed signal [chr", chr_to_plot, "] per replicate and combined")))
}
invisible(dev.off())

# ---------------- media sui bin della maschera (invariato) ----------------
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