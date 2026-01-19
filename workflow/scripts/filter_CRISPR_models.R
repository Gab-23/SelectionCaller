options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args


input_path              <- args$input_path
output_path             <- args$output_path
oncotree_code_DepMap    <- args$oncotree_code_DepMap

models <- read.delim(input_path, sep = ",", header = T)

models_filtered <- models[models$OncotreeCode == oncotree_code_DepMap, "ModelID"]

writeLines(text = models_filtered, con = output_path, sep = "\n")