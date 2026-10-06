options(scipen=999)                                                               # unable scientific notation

# Gene -> interactor table (long format), assembly- and backbone-independent.
#   CORUM: two genes are complex partners if they are subunits of the same complex. Complexes
#          are grouped by ComplexID (PENNE grouped by ComplexName, which is not unique: 2916
#          complexes share 2700 names, so distinct complexes were merged).
#   Interactome INSIDER: undirected PPIs, homodimers removed (as in PENNE).
# A gene is never its own partner. Output columns: gene, partner, source; CORUM subunits
# without any partner appear once with partner = NA.

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_CORUM                <- args$input_path_CORUM
input_path_InteractomeINSIDER   <- args$input_path_InteractomeINSIDER
output_path                     <- args$output_path

load_object <- function(path, name){
  env <- new.env()
  load(path, envir = env)
  if(!exists(name, envir = env, inherits = FALSE)){
    stop(paste0("object '", name, "' not found in ", path))
  }
  get(name, envir = env)
}

# --- CORUM ---------------------------------------------------------------------------------
complexes <- load_object(input_path_CORUM, "human_proteincomplex")
subunits  <- strsplit(as.character(complexes$subunits.Gene.name.), ";")
corum_long <- data.frame(ComplexID = rep(complexes$ComplexID, lengths(subunits)),
                         gene = trimws(unlist(subunits)),                           # PENNE: str_trim
                         stringsAsFactors = FALSE)
corum_long <- unique(corum_long[!is.na(corum_long$gene) & corum_long$gene != "", ])
corum_pairs <- merge(corum_long, corum_long, by = "ComplexID", suffixes = c("", "_partner"))
corum_pairs <- unique(corum_pairs[corum_pairs$gene != corum_pairs$gene_partner, c("gene", "gene_partner")])
colnames(corum_pairs) <- c("gene", "partner")
# Subunits of single-protein complexes have no partner but are still complex subunits
# (they count in complex_density, as in PENNE): keep them with partner = NA.
lonely <- setdiff(unique(corum_long$gene), corum_pairs$gene)
corum_pairs <- rbind(corum_pairs, data.frame(gene = lonely, partner = rep(NA_character_, length(lonely)),
                                             stringsAsFactors = FALSE))
corum_pairs$source <- "CORUM"

# --- Interactome INSIDER --------------------------------------------------------------------
ppis <- load_object(input_path_InteractomeINSIDER, "HS_interfaces_genenames")
# P1/P2 are stored as factors in the RData: convert explicitly (c() on factors returns the
# integer codes in R < 4.1, which would silently turn symbols into numbers)
ppis$P1 <- as.character(ppis$P1)
ppis$P2 <- as.character(ppis$P2)
ppis <- ppis[!is.na(ppis$P1) & !is.na(ppis$P2) & ppis$P1 != "" & ppis$P2 != "", c("P1", "P2")]
ppis <- ppis[ppis$P1 != ppis$P2, ]                                                  # remove homodimers
ppi_pairs <- unique(data.frame(gene    = c(ppis$P1, ppis$P2),                      # undirected
                               partner = c(ppis$P2, ppis$P1),
                               stringsAsFactors = FALSE))
ppi_pairs$source <- "InteractomeINSIDER"

out <- rbind(corum_pairs, ppi_pairs)
out <- out[order(out$source, out$gene, out$partner), ]

write.table(
    x = out,
    file = output_path,
    col.names = T,
    row.names = F,
    quote = F,
    sep = "\t"
)