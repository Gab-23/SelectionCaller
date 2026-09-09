source("workflow/scripts/helpers/import_table.R")

options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path             <- args$input_path
output_path            <- args$output_path

`%>%` <- magrittr::`%>%`

file_names <- readLines(input_path)

epigenome_list <- lapply(X = file_names, FUN = import_table)

average_epigenome <- dplyr::bind_rows(epigenome_list) %>%                          # columns are matched by names, if a column is missing, filled by NA
                     dplyr::mutate(
                        dplyr::across(-binID, ~ tidyr::replace_na(., 0))          # a state column missing for one epigenome means that epigenome has 0bp of
                     ) %>%                                                        # that state (ChromHMM states partition the whole genome), so it must be filled
                                                                                  # with 0 BEFORE averaging, not after: doing it after meant that mean() (with its
                                                                                  # default na.rm=FALSE) silently propagated NA -> the state's average collapsed
                                                                                  # to 0 for every bin as soon as a single epigenome lacked that state column,
                                                                                  # discarding real signal from every other epigenome that did have it.
                     dplyr::group_by(binID) %>%                                   # group by binID
                     dplyr::summarise(
                        dplyr::across(
                            dplyr::everything(), mean),                           # take the average 
                        .groups = 'drop') %>%
                     dplyr::arrange(
                        as.integer(sub(x = binID, "_(.*)", "")),                  # order for chr, bin
                        as.integer(sub(x = binID, "(.*?)_", ""))) %>%
                     dplyr::select(
                        binID,
                        all_of(setdiff(colnames(.), "binID")[order(as.integer(sub(x = setdiff(colnames(.), "binID"), "_(.*)", "")))]) # order for chromatin state
                        ) 

write.table(x = average_epigenome,                                                 # write output
            file = output_path, 
            quote = F, 
            sep = "\t", 
            row.names = F, 
            col.names = T)