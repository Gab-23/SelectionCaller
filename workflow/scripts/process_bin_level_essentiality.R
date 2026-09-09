options(scipen=999)     

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)

input_path              <- args$input_path
output_path             <- args$output_path

`%>%` <- magrittr::`%>%`

bin_level_essentiality <- read.table(input_path, header = F, sep = "\t")
colnames(bin_level_essentiality) <- c("binID", "ess_score")

bin_level_essentiality %>% dplyr::group_by(binID) %>% dplyr::summarise(ess_score = mean(ess_score)) -> aggregated_bin_level_essentiality

write.table(x = aggregated_bin_level_essentiality,
            file = output_path,
            col.names = T,
            row.names = F,
            quote = F,
            sep = "\t")