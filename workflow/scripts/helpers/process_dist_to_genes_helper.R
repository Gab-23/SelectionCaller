process_dist_to_genes <- function(df){
    df$dist <- ifelse(df$dist == -1, NA, df$dist)
    return(df)
}