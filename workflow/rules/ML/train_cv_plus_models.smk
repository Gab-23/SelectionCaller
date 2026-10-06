rule train_cv_plus_models:
    # Second K-fold, usa gli iperparametri/feature consensus fissi (nessuna CV interna).
    #
    # CONFORMAL INTERVALS: abbandonati (2026-09) - la formula CV+ letterale richiederebbe
    # valutare modelli di altri fold su punti che quei fold hanno gia' visto in training
    # (dato che si scora sempre il dataset stesso usato per costruire i fold, mai un punto
    # genuinamente nuovo) - questo rompe l'argomento di exchangeability su cui si basa la
    # garanzia del paper (Barber, Candes, Ramdas & Tibshirani 2021). Sostituiti con un test
    # di significativita' per-bin: z-score robusto (MAD) sui residui out-of-fold, p-value
    # parametrico (coda della normale), correzione BH-FDR - stesso linguaggio statistico gia'
    # usato altrove nella pipeline (get_DARs.smk/edgeR, soglie FDR in COSMIC/GDC).
    #
    # {fdr_threshold} (ex {confidence}) resta un wildcard vero - ogni soglia FDR richiesta fa
    # rigirare l'intero K-fold training/SHAP.
    #
    # {target_col} aggiunto insieme al resto della catena (train_models -> select_best_model)
    # per permettere target diversi da "logFC" senza rifare il DAG a mano.
    #
    # {training_mode} (per_chrom | genome_wide): in genome_wide {chrom} == "all" and the
    # dataset is the full feature matrix (get_training_dataset bypasses the split).
    output:
        models=[
            "results/ML/{dataset_id}/{training_mode}/{chrom}/target/{target_col}/significance/{fdr_threshold}/cv_plus_models/" + str(i) + ".pkl"
            for i in range(OUTER_CV)
        ],
        shap_values="results/ML/{dataset_id}/{training_mode}/{chrom}/target/{target_col}/significance/{fdr_threshold}/shap/shap_values.tsv",
        shap_plots="results/ML/{dataset_id}/{training_mode}/{chrom}/target/{target_col}/significance/{fdr_threshold}/shap/shap_plots.pdf",
        base_value="results/ML/{dataset_id}/{training_mode}/{chrom}/target/{target_col}/significance/{fdr_threshold}/shap/base_value.tsv",
        significance_calls="results/ML/{dataset_id}/{training_mode}/{chrom}/target/{target_col}/significance/{fdr_threshold}/significance_calls.tsv",
        warning="results/ML/{dataset_id}/{training_mode}/{chrom}/target/{target_col}/significance/{fdr_threshold}/logs/cv_plus_warning.log"
    input:
        dataset=get_training_dataset,
        models=input_path_models,
        best_model_name="results/ML/{dataset_id}/{training_mode}/{chrom}/target/{target_col}/best_model/best_model_name.txt",
        best_blocking="results/ML/{dataset_id}/{training_mode}/{chrom}/target/{target_col}/best_model/best_blocking.txt",
        consensus_hyperparams="results/ML/{dataset_id}/{training_mode}/{chrom}/target/{target_col}/best_model/consensus_hyperparams.tsv",
        consensus_features="results/ML/{dataset_id}/{training_mode}/{chrom}/target/{target_col}/best_model/consensus_features.tsv"
    params:
        model_dir="results/ML/{dataset_id}/{training_mode}/{chrom}/target/{target_col}/significance/{fdr_threshold}/cv_plus_models",
        outer_cv=OUTER_CV,
        model_n_jobs=MODEL_N_JOBS_ML,
        bin_size_bp=BIN_SIZE_BP
    log:
        "logs/train_cv_plus_models/{dataset_id}_{training_mode}_{chrom}_{target_col}_{fdr_threshold}.log"
    conda:
        "../../envs/py_ml_env.yml"
    shell:
        """
        rm -f {output.warning}

        python -m workflow.scripts.train_cv_plus_models \
        --input_path_dataset={input.dataset} \
        --input_path_models={input.models} \
        --input_path_best_model_name={input.best_model_name} \
        --input_path_best_blocking={input.best_blocking} \
        --input_path_consensus_hyperparams={input.consensus_hyperparams} \
        --input_path_consensus_features={input.consensus_features} \
        --target_col={wildcards.target_col} \
        --fdr_threshold={wildcards.fdr_threshold} \
        --outer_cv={params.outer_cv} \
        --model_n_jobs={params.model_n_jobs} \
        --bin_size_bp={params.bin_size_bp} \
        --output_dir_models={params.model_dir} \
        --output_path_shap_values={output.shap_values} \
        --output_path_shap_plots={output.shap_plots} \
        --output_path_base_value={output.base_value} \
        --output_path_significance_calls={output.significance_calls} \
        --output_path_warning={output.warning} 2>{log}
        """