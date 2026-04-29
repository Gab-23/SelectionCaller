find_neighbors <- function(data, idx, K = 1, mode = "spatial_analysis"){
    # TODO: add bin number taking into account
    # NOTA: problem because messes up indexing of spdep
	# added mode argument, for smoothing and spatial analysis
	
	# add chr column if not present
	if (!("chr" %in% colnames(data))){
		data$chr <- as.integer(gsub(pattern = "_[0-9]+", replacement = "", x = data$binID))
	}
	
    data$bin_number <- as.integer(gsub(pattern = "[0-9]+_", replacement = "", x = data$binID))
    
    n <- nrow(data)
	if (mode == "spatial_analysis"){neighbors <- c()}
	else if (mode == "smoothing"){neighbors <- c(idx)}    

    for (k in 1:K){

        prev_row <- data[max(0, (idx-k)),]
		curr_row <- data[idx,]
		next_row <- data[(idx+k),]
	
		prev_chr <- prev_row$chr
		curr_chr <- curr_row$chr
		next_chr <- next_row$chr

#		prev_bin_number <- prev_row$bin_number
#       curr_bin_number <- curr_row$bin_number
#       next_bin_number <- next_row$bin_number
	
	if(length(prev_chr) == 0){prev_chr <- "chrNA"; prev_bin_number <- NA}
	if(is.na(next_chr)){next_chr <- "chrNA"; next_bin_number <- NA}
	
	if((idx > 1) & (prev_chr == curr_chr)){
	  neighbors <- c(neighbors, as.integer(idx-k))
	}
	
	if((idx < n) & (next_chr == curr_chr)){
	  neighbors <- c(neighbors, as.integer(idx+k))
	}

    }
	
    return(sort(neighbors)) 

}

