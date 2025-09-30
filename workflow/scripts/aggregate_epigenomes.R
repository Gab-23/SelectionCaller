source("workflow/scripts/import_epigenome.R")

options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path             <- args$input_path
output_path            <- args$output_path

`%>%` <- magrittr::`%>%`

file_names <- readLines(input_path)

epigenome_list <- lapply(X = file_names, FUN = import_epigenome)

average_epigenome <- dplyr::bind_rows(epigenome_list) %>%                         # columns are matched by names, if a column is missing, filled by NA
                     dplyr::group_by(binID) %>%
                     dplyr::summarise(
                        dplyr::across(
                            dplyr::everything(), mean), 
                        .groups = 'drop') %>%
                     dplyr::arrange(
                        as.integer(sub(x = binID, "_(.*)", "")),
                        as.integer(sub(x = binID, "(.*?)_", ""))) %>%
                     dplyr::select(
                        binID,
                        all_of(setdiff(colnames(.), "binID")[order(as.integer(sub(x = setdiff(colnames(.), "binID"), "_(.*)", "")))])
                        ) 

write.table(x = average_epigenome,                                                 # write output
            file = output_path, 
            quote = F, 
            sep = "\t", 
            row.names = F, 
            col.names = T)