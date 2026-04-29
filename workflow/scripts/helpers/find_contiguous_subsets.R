find_contiguous_subsets <- function(chr_subset){
    # helper fuunction for repliseq smoothing
    breaks <- which(diff(as.numeric(chr_subset$pos)) != 1)
    num_breaks <- length(breaks)
    
    out <- c()
    for (break_v in seq(num_breaks)){
        row_idx <- breaks[break_v]
        if (break_v == 1){
            out <- c(out, 1, row_idx, (row_idx+1))
        }
        else if (break_v == num_breaks){
            out <- c(out, row_idx, (row_idx+1), nrow(chr_subset))
        }
        else{
            out <- c(out, row_idx, (row_idx+1))
        }
    }
    out <- sort(out)
    
    init <- 1
    contiguous_subsets <- list()
    for (out_idx in seq(ceiling(length(out) / 2))){
        next_idx <- init + 1
        range <- out[init:next_idx]
        contiguous_subsets[[out_idx]] <- chr_subset[range[1]:range[2],]
        init <- next_idx + 1
    }
    return(contiguous_subsets)
}