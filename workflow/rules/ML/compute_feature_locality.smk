rule compute_feature_locality:
    # Per-feature "locality" score: log((mean |SHAP| + eps) / (n_peaks + 1)) - importance
    # (mean/median |SHAP| across bins) divided by number of local maxima detected along the
    # genomic signal (scipy.signal.find_peaks or similar, with a prominence threshold scaled
    # to that feature's own SHAP signal to stay scale-invariant across features). Positive =
    # broadly/globally important (high importance, few peaks); negative = locally important
    # (concentrated peaks). The +1/+eps smoothing avoids divide-by-zero for a feature with no
    # detectable peaks at all (treated as maximally global) without a separate edge case.
    #
    # {confidence} is a genuine wildcard here too, purely because SHAP values now live under
    # a {confidence}-specific path (train_cv_plus_models writes identical SHAP values under
    # each confidence level it's run for) - feature locality itself has nothing to do with
    # confidence level, this just picks which (redundant, identical) copy of shap_values.tsv
    # to read.
    output:
        "results/ML/{dataset_id}/per_chrom/{chrom}/conformal/{confidence}/locality/feature_locality.tsv"
    input:
        dataset="results/ML/{dataset_id}/per_chrom/{chrom}/dataset.tsv",
        shap_values="results/ML/{dataset_id}/per_chrom/{chrom}/conformal/{confidence}/shap/shap_values.tsv"
    log:
        "logs/compute_feature_locality/{dataset_id}_{chrom}_{confidence}.log"
    conda:
        "../../envs/py_ml_env.yml"
    shell:
        # TODO (user)
        """
        """