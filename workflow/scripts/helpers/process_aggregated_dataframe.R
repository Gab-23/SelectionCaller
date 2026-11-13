process_aggregated_dataframe <- function(df){
    # NOTE (fix, dplyr-only rewrite): sostituisce tibble::rownames_to_column + tidyr::separate,
    # non disponibili in questo env (solo dplyr). I nomi dei modelli reali (MODELS file) non
    # contengono mai underscore, quindi uno strsplit("_") semplice su "Model" produce sempre
    # esattamente 3 pezzi (model_name, blocking, fold), stesso comportamento di tidyr::separate
    # in questo caso specifico.
    blocking_levels <- c("noblocking", "5", "10", "15", "20", "chromosome")
    if (!"Model" %in% colnames(df)) {
        df$Model <- rownames(df)
        rownames(df) <- NULL
        df <- df[, c("Model", setdiff(colnames(df), "Model")), drop = FALSE]
    }

    split_cols <- do.call(rbind, strsplit(df$Model, "_"))
    df$model_name <- split_cols[, 1]
    df$blocking <- factor(split_cols[, 2], levels = blocking_levels)
    df$fold <- as.numeric(split_cols[, 3])
    df$Model <- NULL

    return(df)
}