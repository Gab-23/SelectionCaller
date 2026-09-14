create_ordered_y <- function(data, blocking_col, fold_col) {
    folds <- sort(unique(data[[fold_col]]))
    blocking_levels <- c("noblocking", "5", "10", "15", "20")
    ordered_levels <- tidyr::expand_grid(
      blocking = factor(blocking_levels, levels = blocking_levels),
      fold = folds
    ) %>%
      dplyr::mutate(label = paste0(blocking, ".", fold)) %>%
      dplyr::pull(label)    
    data %>%
      dplyr::mutate(
        y_label = factor(paste0(blocking, ".", fold), levels = ordered_levels)
      ) %>%
      dplyr::pull(y_label)
  }

visualize_models <- function(df_features, df_performances, df_hyperparams){
  # NOTE (fix): df_features/df_performances/df_hyperparams are now expected to already be the
  # output of process_aggregated_dataframe() (called upstream in select_best_model.R) - i.e.
  # model_name/blocking/fold are already real columns (blocking already an ordered factor,
  # fold already numeric), and the original "Model" rowname-derived column no longer exists
  # (process_aggregated_dataframe drops it). The rownames_to_column("Model") -> separate() ->
  # mutate(blocking, fold) steps that used to happen again here were entirely redundant with
  # that upstream processing, and would have broken outright (no "Model" column left to pivot
  # around).

  heatmap_df_features <- df_features %>%
    tidyr::pivot_longer(cols = -c(model_name, blocking, fold), names_to = "Feature", values_to = "Selected")

  heatmap_df_features$y_ordered <- create_ordered_y(heatmap_df_features, "blocking", "fold")

  heatmap_df_performances <- df_performances %>%
    tidyr::pivot_longer(cols = -c(model_name, blocking, fold), names_to = "Metric", values_to = "Value")

  heatmap_df_performances$y_ordered <- create_ordered_y(heatmap_df_performances, "blocking", "fold")

  plot_data <- df_performances %>%
    tidyr::pivot_longer(cols = tidyselect::starts_with("train_"), names_to = "metric", values_to = "value")

  plot_medians <- plot_data %>%
    dplyr::group_by(metric, model_name, blocking) %>%
    dplyr::summarise(median_value = median(value, na.rm = TRUE), .groups = "drop")

  heatmap_features <- ggplot2::ggplot(heatmap_df_features, ggplot2::aes(x = Feature, y = y_ordered, fill = factor(Selected))) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_manual(values = c("0" = "white", "1" = "steelblue"), name = "Selected") +
    ggplot2::facet_wrap(~ model_name, scales = "free_y") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1), panel.grid = ggplot2::element_blank()) +
    ggplot2::labs(title = "Binary Feature Heatmap (Grouped by Model)", x = "Features", y = "Blocking_Fold")

  heatmap_performances <- ggplot2::ggplot(heatmap_df_performances, ggplot2::aes(x = Metric, y = y_ordered, fill = Value)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient2(low = "steelblue", mid = "white", high = "darkred", midpoint = median(heatmap_df_performances$Value, na.rm = TRUE)) +
    ggplot2::facet_wrap(~ model_name, scales = "free_y") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1), panel.grid = ggplot2::element_blank()) +
    ggplot2::labs(title = "Model Performance Heatmap (Grouped by Model)", x = "Metric", y = "Blocking_Fold", fill = "Score")
  
  median_boxplots <- ggplot2::ggplot(plot_data, ggplot2::aes(x = blocking, y = value, fill = model_name)) +
    ggplot2::geom_line(data = plot_medians, ggplot2::aes(x = blocking, y = median_value, color = model_name, group = model_name), position = ggplot2::position_dodge(width = 0.8), linewidth = 1, alpha = 0.7) +
    ggplot2::geom_point(data = plot_medians, ggplot2::aes(x = blocking, y = median_value, color = model_name, group = model_name), position = ggplot2::position_dodge(width = 0.8), size = 2) +
    ggplot2::geom_boxplot(ggplot2::aes(group = interaction(blocking, model_name)), position = ggplot2::position_dodge(width = 0.8), width = 0.7, alpha = 0.8) +
    ggplot2::facet_wrap(~metric, scales = "free_y", ncol = 3) +
    ggplot2::labs(x = "Blocking", y = "Metric value", fill = "Model", color = "Model") +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"), axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  # NOTE: hyperparameter names/types differ per model (see MODELS/GET_HYPERPARAM_SPACE.py) - a
  # single heatmap with one shared color scale across mixed numeric/categorical values wouldn't
  # be meaningful, so this renders the raw per-fold hyperparameter values as a plain text table
  # instead, faceted by model via one grob per model_name. If you had a different visualization
  # in mind for this, let me know and I'll change it.
  hyperparams_table <- df_hyperparams %>%
    dplyr::arrange(model_name, blocking, fold) %>%
    (function(df) {
      tables_per_model <- lapply(split(df, df$model_name), function(df_model) {
        df_model <- df_model[, colSums(!is.na(df_model)) > 0, drop = FALSE]
        gridExtra::tableGrob(df_model, rows = NULL)
      })
      gridExtra::arrangeGrob(grobs = tables_per_model, ncol = 1)
    })()

  to.return <- list(heatmap_features  = heatmap_features, 
                    heatmap_performances = heatmap_performances,
                    median_boxplots = median_boxplots,
                    hyperparams_table = hyperparams_table)

  return(to.return)
}