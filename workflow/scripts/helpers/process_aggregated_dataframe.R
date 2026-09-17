process_aggregated_dataframe <- function(df){
    # NOTE (fix): now tolerates a dataframe that already has a "Model" column (rather than
    # requiring it live in rownames) - needed for the genome-wide pooled performances case in
    # select_best_model.R, where rownames alone are NOT safe to combine across chromosomes:
    # confirmed directly that do.call(rbind, ...) silently mangles duplicate rownames (e.g.
    # two chromosomes both producing "XGBRegressor_5_0") into something like
    # "XGBRegressor_5_01", corrupting the fold number once re-parsed. When pooling across
    # chromosomes, the caller now builds an explicit, disambiguated "Model" column (plus a
    # "chrom" column) BEFORE combining, bypassing the rowname collision entirely.
    if (!"Model" %in% colnames(df)) {
        df <- tibble::rownames_to_column(df, "Model")
    }
    df <- df %>%
        tidyr::separate(Model, sep = "_", into = c("model_name", "blocking", "fold"), remove = F) %>%
        dplyr::mutate(blocking = factor(blocking, levels = c("noblocking", "5", "10", "15", "20")), fold = as.numeric(fold)) %>%
        dplyr::select(-Model)    
    return(df)
}