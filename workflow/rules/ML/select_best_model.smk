rule select_best_model:
    # {chrom} wildcard, but the (model, blocking) choice itself is pooled genome-wide: every
    # chromosome's job independently reads performance files from ALL chromosomes (see
    # params.results_dir) to pick the same winning (model, blocking) - avoiding a per-
    # chromosome choice that could confound "difference due to the chromosome" with
    # "difference due to model architecture" when comparing feature importance across
    # chromosomes later (see conversation). This recomputes the same genome-wide answer once
    # per chromosome (cheap - reading small performance .tsv files, not retraining anything),
    # trading a little redundant computation for a single, simpler rule instead of splitting
    # genome-wide selection and per-chromosome consensus into two rules.
    # Once the winning (model, blocking) is known, hyperparameter/feature consensus is
    # computed from ONLY this chromosome's own folds for that combination - this is where the
    # chromosome-specific signal (the rank-importance correlation you've already validated)
    # needs to live.
    output:
        best_model_name="results/ML/{dataset_id}/per_chrom/{chrom}/best_model/best_model_name.txt",
        best_blocking="results/ML/{dataset_id}/per_chrom/{chrom}/best_model/best_blocking.txt",
        consensus_hyperparams="results/ML/{dataset_id}/per_chrom/{chrom}/best_model/consensus_hyperparams.tsv",
        consensus_features="results/ML/{dataset_id}/per_chrom/{chrom}/best_model/consensus_features.tsv",
        path_plot_selfeats="results/ML/{dataset_id}/per_chrom/{chrom}/plots/selfeats_heatmap.pdf",
        path_plot_performance="results/ML/{dataset_id}/per_chrom/{chrom}/plots/performance_heatmap.pdf",
        path_median_boxplots="results/ML/{dataset_id}/per_chrom/{chrom}/plots/median_boxplots.pdf",
        path_plot_hyperparams="results/ML/{dataset_id}/per_chrom/{chrom}/plots/hyperparams_table.pdf"
    input:
        # every chromosome's train_models must be done - not just this one - since the
        # (model, blocking) choice pools performance evidence genome-wide.
        checks=get_all_chrom_train_model_checks
    params:
        results_dir="results/ML/{dataset_id}",
        hyperparams_dir="results/ML/{dataset_id}/per_chrom/{chrom}/hyperparams",
        selfeats_dir="results/ML/{dataset_id}/per_chrom/{chrom}/selfeats"
    log:
        "logs/select_best_model/{dataset_id}_{chrom}.log"
    conda:
        "../../envs/r_plotting.yml"
    shell:
        """
        Rscript workflow/scripts/select_best_model.R \
        --results_dir={params.results_dir} \
        --input_dir_hyperparams={params.hyperparams_dir} \
        --input_dir_selfeats={params.selfeats_dir} \
        --output_path_best_model_name={output.best_model_name} \
        --output_path_best_blocking={output.best_blocking} \
        --output_path_consensus_hyperparams={output.consensus_hyperparams} \
        --output_path_consensus_features={output.consensus_features} \
        --output_path_plot_selfeats={output.path_plot_selfeats} \
        --output_path_plot_performance={output.path_plot_performance} \
        --output_path_median_boxplots={output.path_median_boxplots} \
        --output_path_plot_hyperparams={output.path_plot_hyperparams} 2>{log}
        """