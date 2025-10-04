source("workflow/scripts/import_table.R")

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

output_path          <- args$output_path
observed_output      <- args$observed_output
expected_output      <- args$expected_output

print(observed_output)
print(expected_output)

observed_output <- import_table(observed_output)
expected_output <- import_table(expected_output)

print(observed_output)
print(expected_output)