options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_paralogs       <-  args$input_path_paralogs
input_path_backbone       <-  args$input_path_backbone
output_path_paralogs      <-  args$output_path_paralogs

load(input_path_paralogs)

backbone_annotated <- read.table(input_path_backbone, sep = "\t", header = F)
colnames(backbone_annotated) <- c("chr", "start", "end", "binID", "genes")

features.df <- c() 
for(binID in backbone_annotated$binID){
	genes <- strsplit(backbone_annotated[backbone_annotated$binID == binID,"genes"],",")[[1]]
	if(length(genes) == 0){
		row <- c(binID,rep(NA,5))
	} else {
		df.paralog <- paralogs[paralogs$symbol %in% genes,]
		df.paralog$n_paralogs <- as.numeric(df.paralog$n_paralogs)
	 	n_paralog_bin <- nrow(df.paralog)
	  if(n_paralog_bin != 0){
	  	df.paralog$n_paralog_trans <- unlist(lapply(df.paralog$paralogs, function(x)length(setdiff(strsplit(x,"::")[[1]],genes)))) # find paralogs on different bins
		mean_n_paralogs <- mean(df.paralog$n_paralogs,na.rm = TRUE)
		total_n_paralogs <- sum(df.paralog$n_paralogs,na.rm = TRUE)
		mean_n_paralogs_trans <- mean(df.paralog$n_paralog_trans,na.rm = TRUE)
		total_n_paralogs_trans <- sum(df.paralog$n_paralog_trans,na.rm = TRUE)
		} else {
		mean_n_paralogs <- NA
		total_n_paralogs <- NA
		mean_n_paralogs_trans <- NA
		total_n_paralogs_trans <- NA
	    }
	  row <- c(binID,n_paralog_bin, mean_n_paralogs, total_n_paralogs, mean_n_paralogs_trans, total_n_paralogs_trans)
	}
	features.df <- rbind(features.df,row)
}
features.df <- as.data.frame(features.df)
colnames(features.df) <- c("binID","n_paralogs_bin", "mean_n_paralogs", "total_n_paralogs", "mean_n_paralogs_trans", "total_n_paralogs_trans")
features.df[,2:6] <- sapply(features.df[,2:6], as.numeric)

write.table(
    x = features.df,
    file = output_path_paralogs,
    col.names = T,
    row.names = F,
    quote = F,
    sep = "\t"
)