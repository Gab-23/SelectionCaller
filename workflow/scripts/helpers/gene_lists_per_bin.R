# Function to find gene groups for a binID
gene.lists.per.binID <- function(binID){
  gene.set <- list("genes" = as.character(strsplit(backbone[backbone$binID == binID, "genes"],"::")[[1]]), #Genes on the bin
                   "partners.trans" = as.character(strsplit(backbone[backbone$binID == binID, "Partners.trans"],"::")[[1]]), #Partners in trans
                   "ppis.trans" = as.character(strsplit(backbone[backbone$binID == binID, "PPIs.trans"],"::")[[1]]), #PPIs in trans
                   "all.int.trans" = union(as.character(strsplit(backbone[backbone$binID == binID, "Partners.trans"],"::")[[1]]),
                                           as.character(strsplit(backbone[backbone$binID == binID, "PPIs.trans"],"::")[[1]])) #All interactions
  ) #For each binID, define a new gene sets
  names(gene.set) <- paste(binID, names(gene.set), sep = "_")
  return(gene.set)
  }