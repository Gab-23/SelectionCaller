rule train_cv_plus_models:
    output:
        models=[
            "results/ML/{dataset_id}/per_chrom/{chrom}/significance/{fdr_threshold}/cv_plus_models/" + str(i) + ".pkl"
            for i in range(OUTER_CV)
        ],
        shap_values="results/ML/{dataset_id}/per_chrom/{chrom}/significance/{fdr_threshold}/shap/shap_values.tsv",
        shap_plots="results/ML/{dataset_id}/per_chrom/{chrom}/significance/{fdr_threshold}/shap/shap_plots.pdf",
        base_value="results/ML/{dataset_id}/per_chrom/{chrom}/significance/{fdr_threshold}/shap/base_value.tsv",
        significance_calls="results/ML/{dataset_id}/per_chrom/{chrom}/significance/{fdr_threshold}/significance_calls.tsv"
    input:
        dataset="results/ML/{dataset_id}/per_chrom/{chrom}/dataset.tsv",
        models=input_path_models,
        best_model_name="results/ML/{dataset_id}/per_chrom/{chrom}/best_model/best_model_name.txt",
        best_blocking="results/ML/{dataset_id}/per_chrom/{chrom}/best_model/best_blocking.txt",
        consensus_hyperparams="results/ML/{dataset_id}/per_chrom/{chrom}/best_model/consensus_hyperparams.tsv",
        consensus_features="results/ML/{dataset_id}/per_chrom/{chrom}/best_model/consensus_features.tsv"
    params:
        model_dir="results/ML/{dataset_id}/per_chrom/{chrom}/significance/{fdr_threshold}/cv_plus_models",
        outer_cv=OUTER_CV,
        model_n_jobs=MODEL_N_JOBS_ML,
        bin_size_bp=BIN_SIZE_BP
    log:
        "logs/train_cv_plus_models/{dataset_id}_{chrom}_{fdr_threshold}.log"
    conda:
        "../../envs/py_ml_env.yml"
    shell:
        """
        python -m workflow.scripts.train_cv_plus_models \
        --input_path_dataset={input.dataset} \
        --input_path_models={input.models} \
        --input_path_best_model_name={input.best_model_name} \
        --input_path_best_blocking={input.best_blocking} \
        --input_path_consensus_hyperparams={input.consensus_hyperparams} \
        --input_path_consensus_features={input.consensus_features} \
        --fdr_threshold={wildcards.fdr_threshold} \
        --outer_cv={params.outer_cv} \
        --model_n_jobs={params.model_n_jobs} \
        --bin_size_bp={params.bin_size_bp} \
        --output_dir_models={params.model_dir} \
        --output_path_shap_values={output.shap_values} \
        --output_path_shap_plots={output.shap_plots} \
        --output_path_base_value={output.base_value} \
        --output_path_significance_calls={output.significance_calls} 2>{log}
        """