source("workflow/scripts/compute_state_length.R")

options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

epigenomes_paths       <- args$epigenomes_paths
backbone_path          <- args$backbone_path
output_path            <- args$output_path

`%>%` <- magrittr::`%>%`

output_path_parts <- paste0(strsplit(x = gsub(x = output_path, 
                                         pattern = "resources/epigenomes_processing_log_", 
                                         replacement = ""
                                         ), 
                                    split = "_"
                                    )[[1]][c(-1,-2)], 
                           collapse = "_"
                           )

epigenomes_paths <- readLines(epigenomes_paths)

output_files <- sapply(X = epigenomes_paths, 
                       FUN = function(x){
                            epigenome_id <- strsplit(x = x, split = "_")[[1]][1]
                            epigenome_output_path <- paste0(epigenome_id, "_", output_path_parts, ".tsv")
                            state_length <- compute_state_length(x, backbone_path, epigenome_output_path)
                            return(epigenome_output_path)
                       }, 
                       simplify = TRUE, 
                       USE.NAMES = FALSE
                      )

writeLines(text = output_files, 
           con = output_path, 
           sep = "\n")