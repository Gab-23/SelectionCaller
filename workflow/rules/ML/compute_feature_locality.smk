rule compute_feature_locality:
    output:
        feature_locality="results/ML/{dataset_id}/per_chrom/{chrom}/target/{target_col}/significance/{fdr_threshold}/locality/{stringency}/{smoothing_size}/feature_locality.tsv"
    input:
        shap_values="results/ML/{dataset_id}/per_chrom/{chrom}/target/{target_col}/significance/{fdr_threshold}/shap/shap_values.tsv"
    params:
        apply_smoothing=APPLY_SMOOTHING,
        importance_aggregation=IMPORTANCE_AGGREGATION
    log:
        "logs/compute_feature_locality/{dataset_id}_{chrom}_{target_col}_{fdr_threshold}_{stringency}_{smoothing_size}.log"
    conda:
        "../../envs/py_ml_env.yml"
    shell:
        """
        python -m workflow.scripts.compute_feature_locality \
        --input_path_shap_values={input.shap_values} \
        --stringency={wildcards.stringency} \
        --apply_smoothing={params.apply_smoothing} \
        --smoothing_size={wildcards.smoothing_size} \
        --importance_aggregation={params.importance_aggregation} \
        --output_path_feature_locality={output.feature_locality} 2>{log}
        """