source("workflow/scripts/helpers/import_table.R")
source("workflow/scripts/helpers/export_table.R")
source("workflow/scripts/helpers/cap_distance_feature_helper.R")

options(scipen=999)                                                               # unable scientific notation

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path                  <- args$input_path
dist_cap                    <- as.numeric(args$dist_cap)
output_path                 <- args$output_path
output_path_plot            <- args$output_path_plot

if (is.na(dist_cap) || dist_cap <= 0) stop("dist_cap must be a positive number of bp, got: ", args$dist_cap)

table <- import_table(input_path)                                                 # binID + one distance column (output of process_dist_to_*)
capped_table <- cap_distance_feature(table, dist_cap)

export_table(capped_table, output_path)                                           # same format as the input

plots <- plot_capping_diagnostics(table, capped_table, dist_cap)
pdf(output_path_plot, width = 14, height = 9)                                     # pages: profiles before, profiles after, per-chromosome summary
print(plots$before)
print(plots$after)
print(plots$summary)
invisible(dev.off())