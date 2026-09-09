source("workflow/scripts/helpers/get_outliers.R")

options(scipen=999)                                                               # unable scientific notation           

args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE)                # read args

input_path_counts                     <- args$input_path_counts
input_path_metadata                   <- args$input_path_metadata
input_path_ampl			              <- args$input_path_ampl
input_path_del                        <- args$input_path_del

mode				                  <- args$mode
degree		                          <- as.integer(args$degree)

output_path_corrected_counts          <- args$output_path_corrected_counts
output_path_comparison_plot           <- args$output_path_comparison_plot

`%>%` <- magrittr::`%>%`

if ((mode == "lm") & (degree != 1)){
  degree <- 1
}

counts <- read.table(input_path_counts, header = T, sep = "\t")
metadata <- read.table(input_path_metadata, header = T, sep = "\t")
ampl <- read.table(input_path_ampl, header = T, sep = "\t")
del <- read.table(input_path_del, header = T, sep = "\t")

merged_table <- Reduce(f = function(x,y){merge(x = x, y = y, by = "binID")}, x = list(counts, ampl, del))
num_cols_merged_table <- ncol(merged_table)

corr_df <- do.call(rbind, 
		   lapply(X = seq_along(1:ncol(merged_table))[-c(1,num_cols_merged_table-1,num_cols_merged_table)], 
			  FUN = function(idx){
				  sample_name <- colnames(merged_table)[idx]
				  col <- merged_table[,idx]
				 
				  pearson_ampl  <- cor(col, merged_table$ampl.frequency, method = "pearson", use = "complete.obs")
				  pearson_del <- cor(col, merged_table$del.frequency, method = "pearson", use = "complete.obs")
				  spearman_ampl <- cor(col, merged_table$ampl.frequency, method = "spearman", use = "complete.obs")
				  spearman_del <- cor(col, merged_table$del.frequency, method = "spearman", use = "complete.obs")
				 
				  return(data.frame(pearson_ampl = pearson_ampl,
								    pearson_del = pearson_del,
				 				    spearman_ampl = spearman_ampl,
				 				    spearman_del = spearman_del,
				 				    row.names = sample_name))
				 
				  }
				 )
				)


corr_df$SRA_id <- rownames(corr_df)
rownames(corr_df) <- NULL

metadata_w_corr <- merge(metadata, corr_df, by = "SRA_id")

df_long <- metadata_w_corr %>%
			tidyr::pivot_longer(
			  cols = starts_with("pearson") | starts_with("spearman"),
			  names_to = "correlation_type",
			  values_to = "value"
			)

df_long$condition <- factor(df_long$condition, levels = c("NAT", "tumor"))

nat_outliers <- df_long %>%
		  dplyr::filter(condition == "NAT") %>%
		  dplyr::group_by(correlation_type) %>%
	          dplyr::mutate(is_outlier = get_outliers(value)) %>%
		  dplyr::ungroup() %>%
		  dplyr::filter(is_outlier) %>%
	          dplyr::select(donor_id, correlation_type) %>%
		  dplyr::distinct()

df_to_label <- df_long %>% dplyr::inner_join(nat_outliers, by = c("donor_id", "correlation_type"))

p <- ggplot2::ggplot(df_long, ggplot2::aes(x = condition, y = value, fill = condition)) +
	ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA) +
	ggplot2::geom_line(
	  ggplot2::aes(group = interaction(donor_id, correlation_type)),
		       color = "gray50",
		       alpha = 0.4
		       ) +
	ggplot2::geom_jitter(width = 0.2, alpha = 0.5) +
	ggplot2::facet_wrap(~ correlation_type, scales = "free_y") +
	ggplot2::geom_text(data = df_to_label,
		           ggplot2::aes(label = donor_id),
		           color = "red",
		           vjust = -0.5,
		           size = 3
		) +
	ggplot2::theme_bw() +
	ggplot2::labs(title = "Correlation distributions with paired connections and NAT outliers labeled",
		      x = "Condition",
		      y = "Correlation value"
		)

###

tumor_samples <- metadata[metadata$condition == "tumor", "SRA_id"]
NAT_samples <- metadata[metadata$condition == "NAT", "SRA_id"]

corrected_tumor_counts <- do.call(cbind, 
				  lapply(X = tumor_samples, 
					 FUN = function(sample_name){
						 if (mode == "lm") {
						   formula <- formula <- as.formula(paste(sample_name, "~ ampl.frequency + del.frequency"))
						   model <- summary(lm(formula = formula, data = merged_table))  
						 } else if (mode == "poly") {
						   formula <- as.formula(paste0(sample_name, 
							  		        " ~ ",
							  		        "poly(ampl.frequency, ", degree, ", raw = TRUE) + ", 
							  		        "poly(del.frequency, ", degree, ", raw = TRUE)")
							  		        )
						   model <- summary(lm(formula = formula, data = merged_table))
						 }
						to.return <- data.frame(X = model$residuals)
						colnames(to.return) <- sample_name
						return(to.return)
							  }
							)
						)

corrected_tumor_counts$binID <- merged_table$binID

reassembled_df <- merge(x = merged_table[,c("binID", NAT_samples)], 
			y = corrected_tumor_counts, by = "binID")
						
reassembled_df[,-1][reassembled_df[,-1] < 0] <- 0

ggplot2::ggsave(output_path_comparison_plot, p, width = 7, height = 7)

write.table(x = reassembled_df,
			file = output_path_corrected_counts,
			col.names = T,
			row.names = F,
			sep = "\t", 
			quote = F)
