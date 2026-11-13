source("workflow/scripts/helpers/fraction_plot.R")
source("workflow/scripts/helpers/feature_boxplots.R")
source("workflow/scripts/helpers/boxplot_0s_non0s.R")
source("workflow/scripts/helpers/plot_heatmap_corr.R")


options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_DARs                               <- args$input_path_DARs
input_path_centromeres                        <- args$input_path_centromeres
input_path_telomeres                          <- args$input_path_telomeres
input_path_coding_mutations                   <- args$input_path_coding_mutations
input_path_coding_mutation_counts             <- args$input_path_coding_mutation_counts   # mut_n_events / mut_n_samples / mut_n_sites
input_path_CDS                                <- args$input_path_CDS                      # cds_bp
input_path_coding_sites                       <- args$input_path_coding_sites             # sites_syn / sites_mis / sites_trunc / sites_other_coding
input_path_noncoding_mutations                <- args$input_path_noncoding_mutations
input_path_dist_to_essential                  <- args$input_path_dist_to_essential
input_path_bin_level_essentiality             <- args$input_path_bin_level_essentiality   # NOTE (30/10/2026): ess_score (essenzialita' DepMap per bin)
input_path_epigenomes                         <- args$input_path_epigenomes
input_path_num_fragile_sites                  <- args$input_path_num_fragile_sites
input_path_dist_fragile_sites                 <- args$input_path_dist_fragile_sites
input_path_GC_content                         <- args$input_path_GC_content
input_path_ampl_frequency                     <- args$input_path_ampl_frequency
input_path_del_frequency                      <- args$input_path_del_frequency
input_path_num_overlaps_protein_coding        <- args$input_path_num_overlaps_protein_coding
input_path_num_overlaps_lncRNA                <- args$input_path_num_overlaps_lncRNA
input_path_num_overlaps_miRNA                 <- args$input_path_num_overlaps_miRNA
input_path_num_overlaps_rRNA                  <- args$input_path_num_overlaps_rRNA
input_path_dist_protein_coding                <- args$input_path_dist_protein_coding
input_path_dist_lncRNA                        <- args$input_path_dist_lncRNA
input_path_dist_miRNA                         <- args$input_path_dist_miRNA
input_path_dist_rRNA                          <- args$input_path_dist_rRNA
input_path_TSG                                <- args$input_path_TSG
input_path_OG                                 <- args$input_path_OG
input_path_repliseq                           <- args$input_path_repliseq
input_path_ohnologs                           <- args$input_path_ohnologs
input_path_paralogs                           <- args$input_path_paralogs
input_path_haploinsufficiency                 <- args$input_path_haploinsufficiency
input_path_hic                                <- args$input_path_hic
input_path_interactors                        <- args$input_path_interactors
input_path_GTEx                               <- args$input_path_GTEx

output_path_dataset                           <- args$output_path_dataset
output_path_corr_plot                         <- args$output_path_corr
output_path_NA_frac_plot                      <- args$output_path_NA_frac_plot
output_path_zero_frac_plot                    <- args$output_path_zero_frac_plot
output_path_boxplots                          <- args$output_path_boxplots
output_path_boxplots_scaled                   <- args$output_path_boxplots_scaled
output_path_boxplots_0s_non0s_muts            <- args$output_path_boxplots_0s_non0s_muts
output_path_boxplots_0s_non0s_chromstates     <- args$output_path_boxplots_0s_non0s_chromstates
#output_path_PCA_plot                         <- args$output_path_PCA_plot
output_path_bool_multicollinearity            <- args$output_path_bool_multicollinearity

FDR_thr                                       <- 0.05

`%>%` <- magrittr::`%>%`

# NOTE (29/10/2026): input esclusi dal config (dataset_excluded_inputs, regola combine_dataset) non vengono passati: il loro
# argomento e' NULL e il file non viene letto. Con tutti gli input presenti il comportamento e' identico a prima.
input_paths <- lapply(X = grep(pattern = "input_path", x = ls(), value = TRUE), FUN = get)
input_paths <- Filter(function(p) !is.null(p) && !isTRUE(p) && nzchar(p), input_paths)

data.frames <- lapply(X = input_paths, 
                      FUN = function(input_path){
                        to.return <- read.delim(file = input_path, header = T, sep = "\t", check.names = F)
                        # NOTE (29/10/2026): righe con binID mancante (es. eventi fuori da ogni bin scritti con binID NA) tolte
                        # PRIMA della fusione: merge() fa corrispondere NA con NA e creava righe spurie (271 a 20 kb nel LUAD).
                        # Nel seno erano gia' tolte dal filtro sul logFC, quindi il dataset del seno non cambia.
                        n_na <- sum(is.na(to.return$binID) | to.return$binID == "")
                        if (n_na > 0) message(sprintf("%s: %d righe con binID mancante tolte", basename(input_path), n_na))
                        to.return <- to.return[!(is.na(to.return$binID) | to.return$binID == ""), , drop = FALSE]
                        to.return <- unique(to.return) 
                        return(to.return)
                      })

dataset <- Reduce(x = data.frames, f = function(x,y){merge(x, y, by="binID", sort = F, all = T)})
# NOTE (29/10/2026): il filtro sui bin senza logFC si applica solo se l'accessibilita' e' nel dataset (input DARs non escluso).
# ATTENZIONE (decisione aperta): toglie i bin senza logFC per TUTTE le omiche (anche mutazioni e CNA).
if ("logFC" %in% colnames(dataset)) {
    dataset <- dataset[!is.na(dataset$logFC),] # remove missing target variable
} else {
    message("colonna logFC assente (accessibilita' esclusa): nessun filtro sui bin senza logFC")
}

dataset.scaled <- as.data.frame(apply(X = dataset[-1], MARGIN = 2, FUN = scale))
dataset.scaled$binID <- dataset$binID

# bin assenti dai file di mutazioni/CDS/siti = 0 mutazioni / 0 bp codificanti / 0 siti
idx_cols_to_impute <- grep(pattern = "length_(non)?coding|^mut_n_|^cds_bp$|^sites_", x = colnames(dataset))
dataset[idx_cols_to_impute][is.na(dataset[idx_cols_to_impute])] <- 0
### Add also imputation for haploinsufficiency

na_fraction <- data.frame(variable = names(dataset[,-1]), frac = sapply(dataset[,-1], function(x) mean(is.na(x))))
zero_fraction <- data.frame(variable = names(dataset[,-1]),frac = sapply(dataset[,-1], function(x) mean(x == 0, na.rm = TRUE)))

long_data <- tidyr::pivot_longer(
  dataset,
  cols = -binID,
  names_to = "variable",
  values_to = "value"
)

long_data.scaled <- tidyr::pivot_longer(
  dataset.scaled,
  cols = -binID,
  names_to = "variable",
  values_to = "value"
)

cor.mat <- cor(dataset[,-1], use = "pairwise.complete.obs", method = "pearson")
cor.df <- as.data.frame(as.table(cor.mat))
names(cor.df) <- c("Var1", "Var2", "Correlation")

NA_frac_plot <- fraction_plot(na_fraction, y_axis = "NA")
zero_frac_plot <- fraction_plot(zero_fraction, y_axis = "zero")

boxplots <- feature_boxplots(long_data)
boxplots.scaled <- feature_boxplots(long_data.scaled)
# i grafici "0 contro non 0" mostrano la distribuzione del logFC: senza accessibilita' si scrive un grafico segnaposto
# (i file di output della regola devono comunque esistere)
placeholder_plot <- function(msg) ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0, label = msg) + ggplot2::theme_void()
if ("logFC" %in% colnames(dataset)) {
    boxplots_0s_non0s_muts <- boxplots_0s_non0s(dataset, "length_(non)?coding")
    boxplots_0s_chromstates <- boxplots_0s_non0s(dataset, "^[0-9]+_[aA-zZ]+")
} else {
    boxplots_0s_non0s_muts <- placeholder_plot("logFC not available (accessibility excluded from the dataset)")
    boxplots_0s_chromstates <- boxplots_0s_non0s_muts
}
heatmap <- plot_heatmap_corr(cor.df)

cols_to_test <- colnames(dataset)[-c(1,2)]

cor.test.result <- do.call(rbind,
                      lapply(X = cols_to_test, 
                             FUN = function(col){
                              x <- dataset[,col]
                              do.call(rbind,
                                 lapply(X = cols_to_test, 
                                        FUN = function(col2){
                                        if(col != col2){
                                            y <- dataset[,col2]
                                            test_result <- tryCatch(
                                                cor.test(x = x, y = y, method = "pearson"),
                                                error = function(e) list(estimate = NA, p.value = NA)
                                            )
                                            row <- data.frame(col1 = col,
                                                              col2 = col2, 
                                                              estimate = round(test_result$estimate,3), 
                                                              pval = round(test_result$p.value, 3))
                                          return(row)
                                        }
                                  }))
                        }))

cor.test.result$FDR <- p.adjust(p = cor.test.result$pval, method = "BH")

bool_multicollinearity <- any(cor.test.result$FDR < FDR_thr, na.rm = TRUE)

# sort DF
coords <- do.call(rbind, strsplit(dataset$binID, "_"))
chr_num <- suppressWarnings(as.numeric(dplyr::recode(coords[,1], X = "23", Y = "24")))   # cromosoma X (e Y) in fondo, in ordine
dataset <- dataset[order(chr_num, as.numeric(coords[,2])),]

writeLines(text = as.character(bool_multicollinearity), con = output_path_bool_multicollinearity, sep = "\n")
write.table(x = dataset, file = output_path_dataset, sep = "\t", quote = F, col.names = T, row.names = F)

ggplot2::ggsave(output_path_NA_frac_plot, NA_frac_plot, width = 12, height = 7)
ggplot2::ggsave(output_path_zero_frac_plot, zero_frac_plot, width = 12, height = 7)
ggplot2::ggsave(output_path_boxplots, boxplots, width = 12, height = 7)
ggplot2::ggsave(output_path_boxplots_scaled, boxplots.scaled, width = 12, height = 7)
ggplot2::ggsave(output_path_boxplots_0s_non0s_muts, boxplots_0s_non0s_muts, width = 12, height = 7)
ggplot2::ggsave(output_path_boxplots_0s_non0s_chromstates, boxplots_0s_chromstates, width = 12, height = 7)
ggplot2::ggsave(output_path_corr_plot, heatmap, width = 12, height = 12)