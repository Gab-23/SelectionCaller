source("workflow/scripts/helpers/CPM.R")

options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_dir                           <- args$input_path_dir
input_path_samples                       <- args$input_path_samples
output_path_data                         <- args$output_path_data
output_path_plot_distribution            <- args$output_path_plot_distribution
output_path_plot_landscape               <- args$output_path_plot_landscape

`%>%` <- magrittr::`%>%`

files <- list.files(input_path_dir, recursive = T, full.names = T, pattern = "filtered")
samples <- read.table(file = input_path_samples, header = T, sep = "\t")

sample_types <- c("ES", "LS", "G1")

binIDs_list <- list()
CPM_list <- list()

for (idx in seq(length(sample_types))){
    
    sample_type <- sample_types[[idx]]

    pattern <- samples[samples$sample_type == sample_type, "SRA_id"]
    library_size <- samples[samples$sample_type == sample_type, "library_size"]

    file <- grep(pattern = pattern, x = files, value = T)
    
    coverage <- read.table(file, header = T, sep = "\t")
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

#print(all(ES.f$binID == common_binIDs))
#print(all(LS.f$binID == common_binIDs))
#print(all(G1.f$binID == common_binIDs))

ES_G1_norm_coverage <- ES.f$CPM / G1.f$CPM
LS_G1_norm_coverage <- LS.f$CPM / G1.f$CPM

EL_coverage <- ES_G1_norm_coverage / LS_G1_norm_coverage

log2_EL_coverage <- log2(EL_coverage)

replication_timing <- data.frame(binID = common_binIDs, 
                                 replication.timing = log2_EL_coverage)

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
  ggplot2::aes(x = pos, y = replication.timing)
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
    y = "Replication timing (log2)",
    title = "Replication timing across chromosomes"
  ) +
  ggplot2::theme_bw()

ggplot2::ggsave(
  filename = output_path_plot_distribution,
  plot = hist,
  width = 8,
  height = 6,
  dpi = 300
)

ggplot2::ggsave(
  filename = output_path_plot_landscape,
  plot = landscape,
  width = 12,
  height = 8,
  units = "in"
)

write.table(x = replication_timing,                                                  # write output
            file = output_path_data, 
            quote = F, 
            sep = "\t", 
            row.names = F, 
            col.names = T)