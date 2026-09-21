options(scipen=999)

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

input_dir_scaling_factors      <- args$input_dir_scaling_factors
input_path_metadata            <- args$input_path_metadata
input_path_samples_to_exclude  <- args$input_path_samples_to_exclude
scaling_pattern                 <- args$scaling_pattern
output_path_plot                <- args$output_path_plot

source("workflow/scripts/helpers/plot_CNVkit_lines.R")

`%>%` <- magrittr::`%>%`

# NOTE: a differenza di get_counts_atacseq.R, qui NON leggiamo counts ne' applichiamo il
# filtro scaling_factors$binID %in% counts$binID - quel filtro serviva solo ad allineare
# scaling_factors a counts per l'analisi differenziale a valle, ma questo grafico usa
# esclusivamente scaling_factors + metadata, quindi e' omesso deliberatamente (non e'
# una dimenticanza) per limitarsi a cio' che serve davvero a rigenerare il grafico.

files_scaling_factors <- list.files(input_dir_scaling_factors, recursive = TRUE, full.names = TRUE, pattern = scaling_pattern)

metadata <- read.table(input_path_metadata, header = T, sep = "\t")
samples_to_exclude <- readLines(input_path_samples_to_exclude)

sample_names_scaling_factors <- unlist(lapply(X = files_scaling_factors, FUN = function(x){basename(dirname(x))}))

scaling_factors_list <- lapply(X = seq(length(files_scaling_factors)),
                   FUN = function(idx){
                    file_scaling_factors <- files_scaling_factors[idx]
                    sample_name_scaling_factors <- sample_names_scaling_factors[idx]
                    tab_scaling_factors <- read.table(file_scaling_factors, header = F, sep = "\t")
                    colnames(tab_scaling_factors) <- c("chr", "start", "end", "binID", "segmentID", sample_name_scaling_factors)
                    to.return <- tab_scaling_factors[,c("binID", sample_name_scaling_factors)]
                    to.return <- to.return %>%
                                    dplyr::group_by(binID) %>%
                                    dplyr::summarise(!!sample_name_scaling_factors := mean(.data[[sample_name_scaling_factors]]))
                    return(to.return)
                    })

scaling_factors <- Reduce(f = function(x, y){merge(x,y, by = "binID")}, x = scaling_factors_list)
scaling_factors <- scaling_factors[,!(colnames(scaling_factors) %in% samples_to_exclude)]

CNVkit_lineplot_mean_only <- plot_CNVkit_lines(scaling_factors, metadata)

ggplot2::ggsave(output_path_plot, plot = CNVkit_lineplot_mean_only, width = 14, height = 4)