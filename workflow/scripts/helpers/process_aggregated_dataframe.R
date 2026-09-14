process_aggregated_dataframe <- function(df){
    df <- df %>%
        tibble::rownames_to_column("Model") %>%
        tidyr::separate(Model, sep = "_", into = c("model_name", "blocking", "fold"), remove = F) %>%
        dplyr::mutate(blocking = factor(blocking, levels = c("noblocking", "5", "10", "15", "20")), fold = as.numeric(fold)) %>%
        dplyr::select(-Model)    
    return(df)
}