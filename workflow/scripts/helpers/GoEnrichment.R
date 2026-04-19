go.enrichment <- function (allGenes, deGenes, gotype = c("BP", "MF", "CC"), cutoff = 0.1, 
                           org = "hs", minCatSize = 10, maxCatSize = 1000, maxResults = 100, 
                           mergeOntologies = F, ID = "SYMBOL") 
{
  require("topGO")
  require("AnnotationDbi")
  gotype <- gsub("GO:", "", gotype, fixed = T)
  gotype <- match.arg(gotype, c("CC", "MF", "BP"), several.ok = T)
  if (substr(org, 0, 2) == "hg") 
    org <- "hs"
  if (substr(org, 0, 2) == "mm") 
    org <- "mm"
  org <- match.arg(org, c("hs", "mm"))
  if (org == "hs") {
    require("org.Hs.eg.db")
  } else {
    require("org.Mm.eg.db")
  }
  if (all(is.na(allGenes))) {
    message("Using all genes as a background.")
    if (org == "hs") {
      data("hg38.refseq")
    } else {
      data("mm10.refseq")
    }
    allGenes <- unique(as.character(gtf$gene_id))
  }
  genes <- as.numeric(allGenes %in% deGenes) + 1
  names(genes) <- allGenes
  ll <- list()
  for (ca in gotype) {
    data <- new("topGOdata", 
                ontology = ca, 
                allGenes = genes, 
                geneSel = function(x) {return(x > 1)}, 
                annot = annFUN.org, 
                mapping = ifelse(org == "hs", "org.Hs.eg.db", "org.Mm.eg.db"), 
                ID = ID, 
                nodeSize = minCatSize)
    resultFisher <- runTest(data, algorithm = "classic", statistic = "fisher")
    df <- GenTable(data, classicFisher = resultFisher, orderBy = "classicFisher", topNodes = length(score(resultFisher)))
    df <- df[df$Annotated < maxCatSize, ]
    df$classicFisher <- gsub(x = df$classicFisher, pattern = "< ", replacement = "")
    df$FDR <- p.adjust(df$classicFisher, method = "fdr")
    df <- df[df$FDR <= cutoff, ]
    if (nrow(df) > 0) {
      if (org == "hs") {
        df$genes <- unlist(lapply(X = AnnotationDbi::mget(as.character(df$GO.ID), org.Hs.egGO2ALLEGS), 
                                  FUN = function(x) {
                                          tmp <- unique(as.character(unlist(AnnotationDbi::mget(as.character(x), org.Hs.egSYMBOL))))
                                          return(paste(tmp[which(tmp %in% deGenes)], sep = ", ", collapse = ", "))
                                        }
                                      )
                                    )
      } else {
        df$genes <- unlist(lapply(AnnotationDbi::mget(as.character(df$GO.ID), 
                                                      org.Mm.egGO2ALLEGS), FUN = function(x) {
                                                        tmp <- unique(as.character(unlist(AnnotationDbi::mget(as.character(x), org.Mm.egSYMBOL))))
                                                        return(paste(tmp[which(tmp %in% deGenes)], sep = ", ", collapse = ", "))
                                                      }))
      }
      df <- df[order(df$classicFisher), ]
      if (nrow(df) > maxResults) 
        df <- df[1:maxResults, ]
      if (length(gotype) == 1) 
        return(df)
      df$Enrichment <- df$Significant/df$Expected
      df <- df[, c("GO.ID", "Term", "Annotated", "Significant", "Enrichment", "classicFisher", "FDR", "genes")]
      names(df)[6] <- "Fisher.pvalue"
      class(df) <- c("data.frame")
      ll[[ca]] <- df
    }
  }
  if (length(ll) == 0) {
    message("Nothing found!")
  } else {
    if (length(ll) == 1) 
      return(ll)
    if (!mergeOntologies) 
      return(ll)
    df <- ll[[1]]
    ll[[1]] <- NULL
    for (i in 1:length(ll)) df <- rbind(df, ll[[i]])
    return(df)
  }
}