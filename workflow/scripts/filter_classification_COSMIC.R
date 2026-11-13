options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args


input_path                       <- args$input_path
output_path_phenotype_id         <- args$output_path_phenotype_id
COSMIC_parameters                <- args$COSMIC_parameters

# NOTE (29/10/2026): filtro facoltativo sul SOTTOTIPO istologico (HISTOLOGY_SUBTYPE_1).
# Formati accettati per COSMIC_parameters:
#   sito_istologia                 (come prima: es. breast_carcinoma -> risultato IDENTICO alla versione precedente)
#   sito_istologia_sottotipo       (es. lung_carcinoma_adenocarcinoma -> LUAD)
#   sito__istologia[__sottotipo]   (doppio "_" come separatore, per valori COSMIC che contengono "_",
#                                   es. lung__carcinoma__squamous_cell_carcinoma, skin__malignant_melanoma)
# Con il singolo "_" i valori non possono contenere "_": un parametro come skin_malignant_melanoma verrebbe letto come
# sito skin, istologia malignant, sottotipo melanoma. Per evitare risultati vuoti senza accorgersene, se nessun fenotipo
# corrisponde lo script si ferma e stampa le combinazioni disponibili per quel sito.
if (grepl("__", COSMIC_parameters, fixed = TRUE)) {
    COSMIC_parameters_parts <- strsplit(COSMIC_parameters, "__", fixed = TRUE)[[1]]
} else {
    COSMIC_parameters_parts <- strsplit(COSMIC_parameters, "_", fixed = TRUE)[[1]]
}
if (!(length(COSMIC_parameters_parts) %in% c(2, 3))) {
    stop("COSMIC_parameters deve avere 2 o 3 campi (sito, istologia, [sottotipo]); per valori con '_' usare '__' come separatore. Trovato: ",
         COSMIC_parameters)
}
site <- COSMIC_parameters_parts[1]
histology <- COSMIC_parameters_parts[2]
subtype <- if (length(COSMIC_parameters_parts) == 3) COSMIC_parameters_parts[3] else NA

classification <- read.delim(input_path, sep = "\t", header = T)

classification_filtered <- classification[(classification$PRIMARY_SITE == site) & (classification$PRIMARY_HISTOLOGY == histology),]
if (!is.na(subtype)) {
    classification_filtered <- classification_filtered[classification_filtered$HISTOLOGY_SUBTYPE_1 == subtype,]
}

if (nrow(classification_filtered) == 0) {
    avail <- unique(classification[classification$PRIMARY_SITE == site, c("PRIMARY_SITE", "PRIMARY_HISTOLOGY", "HISTOLOGY_SUBTYPE_1")])
    stop("nessun fenotipo COSMIC per sito='", site, "', istologia='", histology, "'",
         if (!is.na(subtype)) paste0(", sottotipo='", subtype, "'") else "",
         ". Combinazioni disponibili per questo sito (prime 30):\n",
         paste(utils::capture.output(print(utils::head(avail, 30), row.names = FALSE)), collapse = "\n"))
}
message(sprintf("fenotipi COSMIC selezionati: %d (sito %s, istologia %s, sottotipo %s)",
                length(unique(classification_filtered$COSMIC_PHENOTYPE_ID)), site, histology, ifelse(is.na(subtype), "tutti", subtype)))

#write.table(x = classification_filtered,                                          # write output
#            file = output_path_filtered, 
#            quote = F, 
#            sep = "\t", 
#            row.names = F, 
#            col.names = T)

writeLines(text = unique(classification_filtered$COSMIC_PHENOTYPE_ID), sep = "\n", con = output_path_phenotype_id)