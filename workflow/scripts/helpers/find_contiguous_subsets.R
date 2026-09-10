find_contiguous_subsets <- function(chr_subset){
    # helper function for repliseq smoothing
    breaks <- which(diff(as.numeric(chr_subset$pos)) != 1)
    num_breaks <- length(breaks)

    if (num_breaks == 0){                                                        # fix: explicit early return - the whole subset is one contiguous block.
        return(list(chr_subset))                                                 # avoids the seq(0) gotcha (seq(0) == c(1,0), NOT an empty sequence)
    }

    out <- c()
    for (break_v in seq_len(num_breaks)){                                        # fix: seq_len() instead of seq(), never iterates when num_breaks==0
        row_idx <- breaks[break_v]
        if (break_v == 1){                                                       # fix: independent "if"s (not if/else if) so both the opening and closing
            out <- c(out, 1)                                                     # boundary get added even when a single break is simultaneously the first
        }                                                                        # AND the last one (num_breaks == 1) - previously the elseif chain meant
        out <- c(out, row_idx, (row_idx+1))                                      # only ONE of the two boundaries was ever added in that case, leaving the
        if (break_v == num_breaks){                                              # `out` vector one element short and causing an out-of-bounds NA access
            out <- c(out, nrow(chr_subset))                                      # ("NA/NaN argument") a few lines below - confirmed with a toy test.
        }
    }
    out <- sort(out)

    init <- 1
    contiguous_subsets <- list()
    for (out_idx in seq_len(ceiling(length(out) / 2))){                          # fix: seq_len() here too
        next_idx <- init + 1
        range <- out[init:next_idx]
        contiguous_subsets[[out_idx]] <- chr_subset[range[1]:range[2],, drop = FALSE]
        init <- next_idx + 1
    }
    return(contiguous_subsets)
}
