# Function to calculate mean/median expression values
expression.scores <- function(genes){
  df <- exp.data[exp.data$Gene.name %in% genes,]
  if(nrow(df) == 0){ #In the case of no partner, PPIs or no covarage by expression data
    res <- data.frame("Cohorts" = colnames(df)[-c(1,2)],
                      "Mean" = rep(NA, length(colnames(df)[-c(1,2)])),
                      "Median" = rep(NA, length(colnames(df)[-c(1,2)])))}
  else{
    res <- data.frame("Cohorts" = colnames(df)[-c(1,2)],
                      "Mean" = colMeans(df[,-c(1,2)], na.rm = TRUE),
                      "Median" = apply(df[,-c(1,2)],2,median, na.rm = TRUE))}
  return(res)
}