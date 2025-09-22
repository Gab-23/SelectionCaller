define_windows <- function(mode, window_size, step, chrom_sizes_canon){
  
if(mode %in% c("no", "non-overlapping")){  

windows <- GenomicRanges::tileGenome(seqlengths = chrom_sizes_canon, # named numerical vector
                        	           tilewidth = window_size,        # integer
                                     cut.last.tile.in.chrom = T)     # remove last smaller window
 
} else if(mode %in% c("o", "overlapping")){

chrom_spanning_windows <- data.frame(
  seqnames = names(chrom_sizes_canon),
  start = rep(1, length(chrom_sizes_canon)),
  end = chrom_sizes_canon
)

chrom_spanning_windows <- GenomicRanges::GRanges(seqnames = chrom_spanning_windows$seqnames, 
                                                 ranges   = IRanges::IRanges(start = chrom_spanning_windows$start, 
                                                                             end   = chrom_spanning_windows$end))

windows <- GenomicRanges::slidingWindows(x = chrom_spanning_windows, 
                                         width = window_size, 
                                         step = step)

windows <- unlist(windows)


}

windows <- data.frame(windows)
windows$width <- NULL; windows$strand <- NULL
colnames(windows) <- c("chr", "start", "end")

return(windows)

}
