source("workflow/scripts/helpers/find_neighbors.R")

smooth_something <- function(KNN, something, col, type, C, B){

    res <- do.call(
        cbind,
        parallel::mclapply(
        mc.cores = C,
        mc.set.seed = TRUE,
        X = seq(B), 
        FUN = function(b){

            out <- rep(NA, nrow(something))
            
            shuffled_KNN <- KNN[sample.int(length(KNN), length(KNN))]

            for (idx in seq(shuffled_KNN)) {
                shuffled_KNN.idxs <- shuffled_KNN[[idx]]
                unsmoothed <- something[shuffled_KNN.idxs, col] # annot, or any feature
                
                if(type == "char"){
                  best_guess <- names(which.max(table(unsmoothed)))
                } else if (type == "num"){
                  best_guess <- mean(unsmoothed, na.rm = TRUE)
                }

                out[shuffled_KNN.idxs] <- best_guess
            }

            return(out)

            }
        )
    )

    res <- as.data.frame(res)
    colnames(res) <- paste0("smoothed_", seq(B))
    something <- cbind(something, res)

    if(type == "char"){
       stability <- apply(res, 1, function(x){max(table(x)) / length(x)})
    } else if (type == "num"){
       stability <- apply(res, 1, function(x){sd(x, na.rm = TRUE)})
    }

    smoothed_final <- unlist(
                      lapply(X = seq(nrow(something)), 
    						 FUN = function(idx){

    								row <- unlist(something[idx, grep(pattern = "smoothed", x = colnames(something))])	
                                    
                                    if(type == "char"){
                                        best_guess <- names(which.max(table(row)))
                                    } else if (type == "num"){ # handle NAs
                                        best_guess <- mean(row, na.rm = TRUE)
                                    }

 		    						return(best_guess)
 								}
 							)
                        )

    to.return <- list(smoothed = smoothed_final, stability = stability)

    return(to.return)    

}