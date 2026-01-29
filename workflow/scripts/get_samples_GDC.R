options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_clinical                <- args$input_path_clinical
input_path_sample_sheet            <- args$input_path_sample_sheet
output_path_table                  <- args$output_path_table
output_path_barcodes               <- args$output_path_barcodes
output_path_files                  <- args$output_path_files
gender_GDC                         <- args$gender_GDC
cohort_GDC                         <- args$cohort_GDC

clinical <- read.delim(input_path_clinical, header = T)
clinical_unique <- unique(clinical[,c("project.project_id","cases.case_id","cases.submitter_id","demographic.gender")])
clinical.f <- clinical_unique[(clinical_unique$project.project_id == cohort_GDC) & (clinical_unique$demographic.gender == gender_GDC),]

sample_sheet <- read.delim(input_path_sample_sheet, header = T)

merged <- merge(clinical.f, sample_sheet, by.x = "cases.submitter_id", by.y = "Case.ID")

output <- merged[,c("project.project_id", "cases.submitter_id", "cases.case_id", "Sample.ID", "demographic.gender", "File.ID", "File.Name")]

duplicated <- names(table(output$cases.submitter_id)[which(table(output$cases.submitter_id) > 1)])

output.f <- output[!(output$cases.submitter_id %in% duplicated),]

write.table(x = output.f,                                                  # write output
                file = output_path_table, 
                quote = F, 
                sep = "\t", 
                row.names = F, 
                col.names = T)

writeLines(text = output.f$cases.submitter_id, sep = "\n", con = output_path_barcodes)
writeLines(text = output.f$File.ID, sep = "\n", con = output_path_files)