source("workflow/scripts/compute_state_length.R")

options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

epigenomes_paths       <- args$epigenomes_paths
backbone_path          <- args$backbone_path
output_path            <- args$output_path
is_test_call          <- as.logical(args$is_test_call)

`%>%` <- magrittr::`%>%`

if(!(is_test_call)){                                                             # are we calling this script to test, or for real use?

output_path_parts <- paste0(strsplit(x = gsub(x = output_path,                    # split file name into parts
                                         pattern = "resources/epigenomes_processing_log_", 
                                         replacement = ""
                                         ), 
                                    split = "_"
                                    )[[1]][c(-1,-2)], 
                           collapse = "_"
                           )

} else {

output_path_parts <- "dummy_processed"

}

epigenomes_paths <- readLines(epigenomes_paths)                                   # read file with epigenome paths

output_files <- sapply(X = epigenomes_paths,                                      # for each epigenome
                       FUN = function(x){
                            epigenome_id <- strsplit(x = x, split = "_")[[1]][1]  # take the id
                            epigenome_output_path <- paste0(epigenome_id, "_", output_path_parts, ".tsv") # build the output path
                            state_length <- compute_state_length(x, backbone_path, epigenome_output_path) # compute the matrix bin x state
                            return(epigenome_output_path)                                                 # return the path                             
                       }, 
                       simplify = TRUE, 
                       USE.NAMES = FALSE
                      )

writeLines(text = output_files, #write the paths in a file
           con = output_path, 
           sep = "\n")