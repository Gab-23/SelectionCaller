rule aggregate_and_visualize_models:
    output:
        performance=f"results/ML/{dataset_id}/plots/performance_heatmap.pdf",
        selfeats=f"results/ML/{dataset_id}/plots/selfeats_heatmap.pdf",
        median_boxplots=f"results/ML/{dataset_id}/plots/median_boxplots.pdf"
    input:
        # explicit per-(model,blocking) files, so Snakemake actually waits for every
        # train_model job to finish before aggregating - the R script itself only takes
        # directories (it globs by filename pattern internally), so these aren't referenced
        # directly in the shell command, but they are what makes the dependency real.
        performances=[
            f"results/ML/{dataset_id}/performances/perf_file_{model}_{blocking}.tsv"
            for model in MODELS for blocking in BLOCKINGS
        ],
        selfeats=[
            f"results/ML/{dataset_id}/selfeats/selfeats_file_{model}_{blocking}.tsv"
            for model in MODELS for blocking in BLOCKINGS
        ],
        hyperparams=[
            f"results/ML/{dataset_id}/hyperparams/hyperparams_{model}_{blocking}.tsv"
            for model in MODELS for blocking in BLOCKINGS
        ]
    params:
        input_dir_performance=f"results/ML/{dataset_id}/performances",
        input_dir_selfeats=f"results/ML/{dataset_id}/selfeats",
        # NOTE: input_dir_hyperparam is accepted as a CLI argument by
        # aggregate_and_visualize_models.R (and is required, so we still pass it), but it is
        # never actually referenced anywhere in that script's body - dead/unused input,
        # confirmed by inspection. Kept here only because the script requires the argument.
        input_dir_hyperparam=f"results/ML/{dataset_id}/hyperparams"
    log:
        f"logs/aggregate_and_visualize_models/{dataset_id}.log"
    conda:
        "../../envs/visualize_models.yml"
    shell:
        """
        Rscript workflow/scripts/aggregate_and_visualize_models.R \
        --input_dir_performance={params.input_dir_performance} \
        --input_dir_selfeats={params.input_dir_selfeats} \
        --input_dir_hyperparam={params.input_dir_hyperparam} \
        --output_path_performance={output.performance} \
        --output_path_selfeats={output.selfeats} \
        --output_path_median_boxplots={output.median_boxplots} 2>{log}
        """
