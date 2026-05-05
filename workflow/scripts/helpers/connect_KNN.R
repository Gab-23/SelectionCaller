connect_KNN <- function(KNN){
    
    i <- 1
    j <- i + 1
    k <- i
    cluster_name <- i
    cluster <- KNN[[k]]

    KNN.connected <- list()

    while (j < length(KNN)){
        
        while(length(intersect(cluster, KNN[[j]])) > 0){
            cluster <- sort(unique(unlist(KNN[k:j])))
            i <- i + 1
            j <- i + 1
        }

        KNN.connected[[cluster_name]] <-  cluster
        
        k <- j
        i <- k
        j <- i + 1
        cluster_name <- cluster_name + 1
        cluster <- KNN[[k]]
        
    }

    KNN.connected[[cluster_name]] <-  cluster
    
    return(KNN.connected)
}