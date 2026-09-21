rule select_best_model:
    output:
        best_model_name="results/ML/{dataset_id}/per_chrom/{chrom}/target/{target_col}/best_model/best_model_name.txt",
        best_blocking="results/ML/{dataset_id}/per_chrom/{chrom}/target/{target_col}/best_model/best_blocking.txt",
        consensus_hyperparams="results/ML/{dataset_id}/per_chrom/{chrom}/target/{target_col}/best_model/consensus_hyperparams.tsv",
        consensus_features="results/ML/{dataset_id}/per_chrom/{chrom}/target/{target_col}/best_model/consensus_features.tsv",
        path_plot_selfeats="results/ML/{dataset_id}/per_chrom/{chrom}/target/{target_col}/plots/selfeats_heatmap.pdf",
        path_plot_performance="results/ML/{dataset_id}/per_chrom/{chrom}/target/{target_col}/plots/performance_heatmap.pdf",
        path_median_boxplots="results/ML/{dataset_id}/per_chrom/{chrom}/target/{target_col}/plots/median_boxplots.pdf",
        path_plot_hyperparams="results/ML/{dataset_id}/per_chrom/{chrom}/target/{target_col}/plots/hyperparams_table.pdf"
    input:
        checks=get_all_chrom_train_model_checks
    params:
        results_dir="results/ML/{dataset_id}",
        target_col="{target_col}",
        hyperparams_dir="results/ML/{dataset_id}/per_chrom/{chrom}/target/{target_col}/hyperparams",
        selfeats_dir="results/ML/{dataset_id}/per_chrom/{chrom}/target/{target_col}/selfeats"
    log:
        "logs/select_best_model/{dataset_id}_{chrom}_{target_col}.log"
    conda:
        "../../envs/r_plotting.yml"
    shell:
        """
        Rscript workflow/scripts/select_best_model.R \
        --results_dir={params.results_dir} \
        --target_col={params.target_col} \
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