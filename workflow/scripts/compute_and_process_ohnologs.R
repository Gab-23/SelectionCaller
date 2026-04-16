options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_ohnologs       <-  args$input_path_ohnologs
input_path_backbone       <-  args$input_path_backbone
output_path_ohnologs      <-  args$output_path_ohnologs

load(input_path_ohnologs)

backbone_annotated <- read.table(input_path_backbone, sep = "\t", header = F)
colnames(backbone_annotated) <- c("chr", "start", "end", "binID", "genes")

features.df <- c() 
for(binID in backbone_annotated$binID){
	genes <- strsplit(backbone_annotated[backbone_annotated$binID == binID,"genes"],",")[[1]]
	if(length(genes) == 0){
		row <- c(binID,rep(NA,5))
	} else {
	  	ohnolog <- as.data.frame(ohnologs[["mmpaper"]])
	  	df.ohnolog <- ohnolog[ohnolog$Symbol %in% genes,]
	  	n_ohnolog_bin <- nrow(df.ohnolog)
	  if(n_ohnolog_bin != 0){
		df.ohnolog$n_ohnolog.mmpaper_trans <- unlist(lapply(df.ohnolog$ohnologs_mmpaper, function(x)length(setdiff(strsplit(x,"::")[[1]],genes))))
		mean_n_ohnologs.mmpaper <- mean(df.ohnolog$n_ohnolog_mmpaper,na.rm = TRUE)
		total_n_ohnologs.mmpaper <- sum(df.ohnolog$n_ohnolog_mmpaper,na.rm = TRUE)
		mean_n_ohnologs.mmpaper_trans <- mean(df.ohnolog$n_ohnolog.mmpaper_trans,na.rm = TRUE)
		total_n_ohnologs.mmpaper_trans <- sum(df.ohnolog$n_ohnolog.mmpaper_trans,na.rm = TRUE)
	  } else {
		mean_n_ohnologs.mmpaper <- NA
		total_n_ohnologs.mmpaper <- NA
		mean_n_ohnologs.mmpaper_trans <- NA
		total_n_ohnologs.mmpaper_trans <- NA
	  }
	  row <- c(binID,n_ohnolog_bin,mean_n_ohnologs.mmpaper, total_n_ohnologs.mmpaper, mean_n_ohnologs.mmpaper_trans, total_n_ohnologs.mmpaper_trans)
  }
features.df <- rbind(features.df,row)
}
features.df <- as.data.frame(features.df)
colnames(features.df) <- c("binID","n_ohnolog_bin", "mean_n_ohnologs.mmpaper", "total_n_ohnologs.mmpaper", "mean_n_ohnologs.mmpaper_trans", "total_n_ohnologs.mmpaper_trans")
features.df[,2:6] <- sapply(features.df[,2:6], as.numeric)

write.table(
    x = features.df,
    file = output_path_ohnologs,
    col.names = T,
    row.names = F,
    quote = F,
    sep = "\t"
)