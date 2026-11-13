source("workflow/scripts/helpers/plot_landscape.R")

options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_ampl                     <- args$input_path_ampl
input_path_del                      <- args$input_path_del
output_path_ampl                    <- args$output_path_ampl
output_path_del                     <- args$output_path_del
output_path_plot_ampl               <- args$output_path_plot_ampl
output_path_plot_del                <- args$output_path_plot_del
num_patients                        <- as.numeric(args$num_patients)

`%>%` <- magrittr::`%>%`

ampl <- read.delim(input_path_ampl, header = F)
colnames(ampl) <- c("chr", "start", "end", "binID", "ampl.frequency")
ampl$ampl.frequency <- (ampl$ampl.frequency / num_patients)

del <- read.delim(input_path_del, header = F)
colnames(del) <- c("chr", "start", "end", "binID", "del.frequency")
del$del.frequency <- (del$del.frequency / num_patients)

# NOTE (26/10/2026): aggiunto il cromosoma X ai livelli (prima 1:22: l'X diventava NA nel fattore, quindi spariva dal
# grafico; i file di output, che contengono solo binID e frequenza, non erano toccati).
chr_levels <- c(1:22, "X")

ampl.out <- ampl %>%
  dplyr::select(binID, ampl.frequency) %>%
  dplyr::mutate(
    chromosome = sub("_.*", "", binID),
    bin = as.numeric(sub(".*_", "", binID)),
    chromosome = factor(chromosome, levels = chr_levels)
  )

del.out <- del %>%
  dplyr::select(binID, del.frequency) %>%
  dplyr::mutate(
    chromosome = sub("_.*", "", binID),
    bin = as.numeric(sub(".*_", "", binID)),
    chromosome = factor(chromosome, levels = chr_levels)
  )

p_ampl <- plot_landscape(df = ampl.out, 
                         what_to_plot = "ampl.frequency",
                         y_axis = "Amplification frequency",
                         title = "Amplification frequency across genome")

p_del <- plot_landscape(df = del.out, 
                         what_to_plot = "del.frequency",
                         y_axis = "Deletion frequency",
                         title = "Deletion frequency across genome")

write.table(x = ampl.out[,c("binID", "ampl.frequency")],                                                         # write output
            file = output_path_ampl, 
            quote = F, 
            sep = "\t", 
            row.names = F, 
            col.names = T)

write.table(x = del.out[,c("binID", "del.frequency")],                                                          # write output
            file = output_path_del, 
            quote = F, 
            sep = "\t", 
            row.names = F, 
            col.names = T)

ggplot2::ggsave(output_path_plot_del, p_del, width = 14, height = 4)
ggplot2::ggsave(output_path_plot_ampl, p_ampl, width = 14, height = 4)