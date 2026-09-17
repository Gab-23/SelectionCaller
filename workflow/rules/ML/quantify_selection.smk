rule quantify_selection:
    # Per-bin selection call: check whether the observed value falls within the conformal
    # interval (already computed by train_cv_plus_models, at this specific confidence level)
    # - if yes, the dominant SHAP feature(s) explain the prediction; if no, the sign (and,
    # per the earlier discussion, possibly magnitude normalized by interval width) of the
    # residual is the signal instead.
    #
    # {confidence} is a genuine wildcard (multiple levels explorable) - kept cheap
    # deliberately: this rule only reads the specific confidence level's already-computed
    # intervals.tsv, never touches the K models or SHAP computation.
    output:
        "results/ML/{dataset_id}/per_chrom/{chrom}/selection_calls/{confidence}/selection_calls.tsv"
    input:
        dataset="results/ML/{dataset_id}/per_chrom/{chrom}/dataset.tsv",
        conformal_intervals="results/ML/{dataset_id}/per_chrom/{chrom}/conformal/{confidence}/intervals.tsv",
        shap_values="results/ML/{dataset_id}/per_chrom/{chrom}/conformal/{confidence}/shap/shap_values.tsv"
    log:
        "logs/quantify_selection/{dataset_id}_{chrom}_{confidence}.log"
    conda:
        "../../envs/py_ml_env.yml"
    shell:
        # TODO (user)
        """
        """