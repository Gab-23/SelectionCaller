process_dist_to_cancer_genes <- function(df){
    df$gene <- ifelse(df$gene == ".", NA, df$gene)
    df$role <- ifelse(df$role == ".", NA, df$role)
    df$dist <- ifelse(df$dist == -1, NA, df$dist)
    return(df)
}