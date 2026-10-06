create_ordered_y <- function(data, blocking_col, fold_col) {
    folds <- sort(unique(data[[fold_col]]))
    blocking_levels <- c("noblocking", "5", "10", "15", "20")
    # NOTE (fix, dplyr-only rewrite): tidyr::expand_grid(blocking=..., fold=...) fa variare
    # l'ULTIMO argomento piu' velocemente (blocking costante, fold che cicla) - comportamento
    # OPPOSTO a expand.grid() base, che varia il PRIMO argomento piu' velocemente. Invertendo
    # l'ordine degli argomenti passati a expand.grid() si ottiene la stessa identica sequenza
    # (verificato direttamente).
    grid <- expand.grid(fold = folds, blocking = factor(blocking_levels, levels = blocking_levels))
    ordered_levels <- paste0(grid$blocking, ".", grid$fold)

    y_label <- factor(paste0(data[[blocking_col]], ".", data[[fold_col]]), levels = ordered_levels)
    return(y_label)
}

# NOTE (fix, dplyr-only rewrite): sostituisce tidyr::pivot_longer, non disponibile in questo
# env (solo dplyr, non l'intero tidyverse).
pivot_longer_base <- function(df, id_cols, value_cols, names_to, values_to) {
  do.call(rbind, lapply(value_cols, function(col) {
    out <- df[, id_cols, drop = FALSE]
    out[[names_to]] <- col
    out[[values_to]] <- df[[col]]
    out
  }))
}

visualize_models <- function(df_features, df_performances, df_hyperparams){

  id_cols_features <- c("model_name", "blocking", "fold")
  feature_cols <- setdiff(colnames(df_features), id_cols_features)
  heatmap_df_features <- pivot_longer_base(df_features, id_cols_features, feature_cols, "Feature", "Selected")

  heatmap_df_features$y_ordered <- create_ordered_y(heatmap_df_features, "blocking", "fold")

  id_cols_perf <- c("model_name", "blocking", "fold")
  metric_cols <- setdiff(colnames(df_performances), id_cols_perf)
  heatmap_df_performances <- pivot_longer_base(df_performances, id_cols_perf, metric_cols, "Metric", "Value")

  heatmap_df_performances$y_ordered <- create_ordered_y(heatmap_df_performances, "blocking", "fold")

  # NOTE (fix, dplyr-only rewrite): sostituisce tidyselect::starts_with("train_")
  train_metric_cols <- grep("^train_", colnames(df_performances), value = TRUE)
  plot_data <- pivot_longer_base(df_performances, id_cols_perf, train_metric_cols, "metric", "value")

  plot_medians <- plot_data %>%
    dplyr::group_by(metric, model_name, blocking) %>%
    dplyr::summarise(median_value = median(value, na.rm = TRUE), .groups = "drop")

  # NOTE (fix): un singolo outlier estremo (es. adj_r2 con n-p-1 vicino a zero su un fold
  # piccolo) trascina l'asse y dell'intero pannello facet. outlier.shape=NA da solo NON basta
  # - ggplot calcola comunque il range dell'asse dai dati sottostanti (verificato
  # direttamente). Serve rimuovere davvero le righe estreme dal dato passato al boxplot, per
  # metrica. plot_medians resta calcolato sui dati COMPLETI non filtrati - la mediana e' gia'
  # robusta agli outlier di per se'.
  outlier_bounds <- plot_data %>%
    dplyr::group_by(metric) %>%
    dplyr::summarise(
      lower = quantile(value, 0.01, na.rm = TRUE),
      upper = quantile(value, 0.99, na.rm = TRUE),
      .groups = "drop"
    )

  plot_data_trimmed <- plot_data %>%
    dplyr::left_join(outlier_bounds, by = "metric") %>%
    dplyr::mutate(is_outlier = value < lower | value > upper)

  n_clipped <- sum(plot_data_trimmed$is_outlier, na.rm = TRUE)
  if (n_clipped > 0) {
    message(sprintf("visualize_models: %d observations (su %d) excluded from boxplot because outside of the range 1%%-99%% - median is computed on complete data.", n_clipped, nrow(plot_data_trimmed)))
  }

  plot_data_trimmed <- plot_data_trimmed %>%
    dplyr::filter(!is_outlier)

  heatmap_features <- ggplot2::ggplot(heatmap_df_features, ggplot2::aes(x = Feature, y = y_ordered, fill = factor(Selected))) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_manual(values = c("0" = "white", "1" = "steelblue"), name = "Selected") +
    ggplot2::facet_wrap(~ model_name, scales = "free_y") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1), panel.grid = ggplot2::element_blank()) +
    ggplot2::labs(title = "Binary Feature Heatmap (Grouped by Model)", x = "Features", y = "Blocking_Fold")

    # NOTE (fix): stesso problema di root-cause di median_boxplots, mai corretto qui - un
  # singolo outlier estremo (es. adj_r2 catastrofico su un fold piccolo) domina il range di
  # scale_fill_gradient2, facendo collassare a bianco quasi tutte le celle "normali"
  # (verificato direttamente: con un solo outlier su 20 celle, 19 diventano bianco puro).
  # Serve calcolare midpoint/limits SOLO sui valori puliti (esclusi gli outlier), per metrica,
  # poi usare oob=scales::squish per schiacciare i valori fuori range al colore limite invece
  # di lasciare che stirino l'intera scala. Attenzione: i bound vanno calcolati sui dati GIA'
  # ripuliti, non sui quantili dei dati grezzi (un quantile 1%/99% calcolato su pochi punti
  # ancora contaminati dall'outlier resta esso stesso vicino all'outlier - stesso errore
  # commesso e corretto durante lo sviluppo di questo fix).
  perf_bounds <- heatmap_df_performances %>%
    dplyr::group_by(Metric) %>%
    dplyr::mutate(
      q1 = quantile(Value, 0.25, na.rm = TRUE),
      q3 = quantile(Value, 0.75, na.rm = TRUE),
      iqr = q3 - q1,
      lo = q1 - 3 * iqr,
      hi = q3 + 3 * iqr
    ) %>%
    dplyr::ungroup()

  clean_values <- perf_bounds$Value[perf_bounds$Value >= perf_bounds$lo & perf_bounds$Value <= perf_bounds$hi]

  heatmap_performances <- ggplot2::ggplot(heatmap_df_performances, ggplot2::aes(x = Metric, y = y_ordered, fill = Value)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient2(
      low = "steelblue", mid = "white", high = "darkred",
      midpoint = median(clean_values, na.rm = TRUE),
      limits = range(clean_values, na.rm = TRUE),
      oob = scales::squish
    ) +
    ggplot2::facet_wrap(~ model_name, scales = "free_y") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1), panel.grid = ggplot2::element_blank()) +
    ggplot2::labs(title = "Model Performance Heatmap (Grouped by Model)", x = "Metric", y = "Blocking_Fold", fill = "Score")

    median_boxplots <- ggplot2::ggplot(plot_data_trimmed, ggplot2::aes(x = blocking, y = value, fill = model_name)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.4) +
    ggplot2::geom_line(data = plot_medians, ggplot2::aes(x = blocking, y = median_value, color = model_name, group = model_name), position = ggplot2::position_dodge(width = 0.8), linewidth = 1, alpha = 0.7) +
    ggplot2::geom_point(data = plot_medians, ggplot2::aes(x = blocking, y = median_value, color = model_name, group = model_name), position = ggplot2::position_dodge(width = 0.8), size = 2) +
    ggplot2::geom_boxplot(ggplot2::aes(group = interaction(blocking, model_name)), position = ggplot2::position_dodge(width = 0.8), width = 0.7, alpha = 0.8) +
    ggplot2::facet_wrap(~metric, scales = "free_y", ncol = 3) +
    ggplot2::labs(x = "Blocking", y = "Metric value (1%-99% per metrica, outlier esclusi)", fill = "Model", color = "Model") +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"), axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    # NOTE (nuovo): density plot delle performance, poolate su TUTTI i cromosomi E tutti i
    # modelli (non solo questo cromosoma, non distinto per modello), per confrontare
    # direttamente le strategie di blocking - riusa plot_data_trimmed gia' calcolato sopra
    # (stesso trim degli outlier, stesse colonne model_name/blocking/metric/value).
    blocking_density_plot <- ggplot2::ggplot(plot_data_trimmed, ggplot2::aes(x = value, fill = blocking)) +
      ggplot2::geom_density(alpha = 0.3) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.4) +
      ggplot2::facet_wrap(~ metric, scales = "free", ncol = 3) +
      ggplot2::labs(title = "Performance distribution pooled across chromosomes and models, by blocking strategy",
                    x = "Metric value (1%-99%, outliers excluded)", y = "Density", fill = "Blocking") +
      ggplot2::theme_bw() +
      ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"))

  # NOTE (redesign, sostituisce la vecchia tabella di testo gridExtra::tableGrob): separa
  # iperparametri numerici (jitter per fold, scala libera per iperparametro) da categorici
  # (conteggio a barre) - una tabella di solo testo non mostrava variabilita' o consenso tra
  # fold, e uno scatter/boxplot unificato non avrebbe senso con scale/tipi cosi' diversi tra
  # iperparametri di modelli differenti.
  id_cols_hp <- c("model_name", "blocking", "fold")
  hp_cols <- setdiff(colnames(df_hyperparams), id_cols_hp)

  hp_long <- do.call(rbind, lapply(hp_cols, function(col) {
    out <- df_hyperparams[, id_cols_hp, drop = FALSE]
    out$hyperparam <- sub("^model__", "", col)
    out$value_raw <- as.character(df_hyperparams[[col]])
    out
  }))
  hp_long <- hp_long[!is.na(hp_long$value_raw), ]
  hp_long$value_num <- suppressWarnings(as.numeric(hp_long$value_raw))
  hp_long$is_numeric <- !is.na(hp_long$value_num)

  hp_numeric <- hp_long[hp_long$is_numeric, ]
  hp_categorical <- hp_long[!hp_long$is_numeric, ]

  hp_plots <- list()

  if (nrow(hp_numeric) > 0) {
    hp_plots$numeric <- ggplot2::ggplot(hp_numeric, ggplot2::aes(x = blocking, y = value_num, color = model_name)) +
      ggplot2::geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
      ggplot2::facet_wrap(~ model_name + hyperparam, scales = "free_y") +
      ggplot2::theme_bw() +
      ggplot2::labs(title = "Iperparametri numerici (per fold)", x = "Blocking", y = "Valore", color = "Model")
  }

  if (nrow(hp_categorical) > 0) {
    hp_plots$categorical <- ggplot2::ggplot(hp_categorical, ggplot2::aes(x = value_raw, fill = model_name)) +
      ggplot2::geom_bar(position = "dodge") +
      ggplot2::facet_wrap(~ model_name + hyperparam, scales = "free") +
      ggplot2::theme_bw() +
      ggplot2::labs(title = "Iperparametri categorici (conteggio per fold)", x = "Valore", y = "Conteggio", fill = "Model")
  }

  hp_plots_list <- Filter(Negate(is.null), hp_plots)
  hyperparams_table <- gridExtra::arrangeGrob(grobs = hp_plots_list, ncol = 1)

  to.return <- list(heatmap_features  = heatmap_features, 
                    heatmap_performances = heatmap_performances,
                    median_boxplots = median_boxplots,
                    blocking_density_plot = blocking_density_plot,
                    hyperparams_table = hyperparams_table)

  return(to.return)
}