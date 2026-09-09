source("workflow/scripts/helpers/import_table.R")

options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_IntOGen                         <- args$input_path_IntOGen
input_path_protein_coding_transcript       <- args$input_path_protein_coding_transcript
tumor_type_IntOGen                         <- args$tumor_type_IntOGen
cohort                                     <- args$cohort
output_path                                <- args$output_path

`%>%` <- magrittr::`%>%`

protein_coding_transcript <- import_table(input_path_protein_coding_transcript, header = F)
colnames(protein_coding_transcript) <- c("chr","start","end", "TxID", "gene_name", "gene_type","type")
protein_coding_transcript <- protein_coding_transcript[,c("chr","start","end", "TxID")]

IntOGen <- import_table(input_path_IntOGen)
IntOGen <- IntOGen[,c("SYMBOL", "TRANSCRIPT", "COHORT", "CANCER_TYPE", "ROLE", "IS_DRIVER")]

role_map <- c(Act = "OG", LoF = "TSG", ambiguous = "ambiguous")

unknown_roles <- setdiff(unique(IntOGen$ROLE), names(role_map))
if (length(unknown_roles) > 0){
    stop(sprintf("process_IntOGen.R: unexpected ROLE value(s) not present in the mapping: %s",
                 paste(unknown_roles, collapse = ", ")))
}
IntOGen$ROLE <- unname(role_map[IntOGen$ROLE])
IntOGen$IS_DRIVER <- ifelse(IntOGen$IS_DRIVER == "True", TRUE, FALSE)

IntOGen.f <- IntOGen[(IntOGen$CANCER_TYPE == tumor_type_IntOGen) & (IntOGen$COHORT == cohort) & (IntOGen$IS_DRIVER == TRUE),]
IntOGen.f <- IntOGen.f[,c("SYMBOL", "TRANSCRIPT", "ROLE")]
colnames(IntOGen.f) <- c("gene", "TxID", "role")

IntOGen.f_with_coords <- merge(protein_coding_transcript, IntOGen.f, by = "TxID")
IntOGen.f_with_coords <- IntOGen.f_with_coords[,c("chr", "start", "end", "TxID", "gene", "role")]

write.table(x = IntOGen.f_with_coords,                                                 # write output
            file = output_path, 
            quote = F, 
            sep = "\t", 
            row.names = F, 
            col.names = F)