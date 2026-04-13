process_topTags <- function(topTags){

  bin_split <- do.call(rbind, strsplit(topTags$binID, "_"))
  topTags$chr   <- bin_split[,1]
  topTags$start <- as.numeric(bin_split[,2])

  topTags$chr_num <- gsub("chr", "", topTags$chr)
  topTags$chr_num[topTags$chr_num == "X"] <- 23
  topTags$chr_num[topTags$chr_num == "Y"] <- 24
  topTags$chr_num <- as.numeric(topTags$chr_num)

  topTags <- topTags[order(topTags$chr_num, topTags$start), ]
  topTags$significant <- topTags$FDR < 0.05
  
  return(topTags)

}
