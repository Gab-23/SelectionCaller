connect_KNN_igraph <- function(KNN) {

  null_idxs <- which(unlist(lapply(KNN, is.null)))
  KNN[null_idxs] <- null_idxs

  edges <- c()
  for (i in seq_along(KNN)) {
    for (j in seq_along(KNN)) {
      bool_overlap <- (length(intersect(KNN[[i]], KNN[[j]])) > 0)
      if (i < j && bool_overlap) {
        edges <- c(edges, i, j) # build adj list if KNN[[i]] and KNN[[j]] share overlap
      }
    }
  }

  g <- igraph::graph(edges = edges, directed = FALSE)
  comps <- igraph::components(g)

  clusters_idxs <- split(x = seq_along(KNN), f = comps$membership)

  KNN.connected <- lapply(clusters_idxs, function(idx){sort(unique(unlist(KNN[idx])))})

  return(KNN.connected)
}