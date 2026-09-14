rule select_best_model:
    output:
        best_model_name="results/ML/{dataset_id}/per_chrom/{chrom}/best_model/best_model_name.txt",
        consensus_hyperparams="results/ML/{dataset_id}/per_chrom/{chrom}/best_model/consensus_hyperparams.tsv",
        consensus_features="results/ML/{dataset_id}/per_chrom/{chrom}/best_model/consensus_features.tsv",
        path_plot_selfeats="results/ML/{dataset_id}/per_chrom/{chrom}/plots/selfeats_heatmap.pdf",
        path_plot_performance="results/ML/{dataset_id}/per_chrom/{chrom}/plots/performance_heatmap.pdf",
        path_median_boxplots="results/ML/{dataset_id}/per_chrom/{chrom}/plots/median_boxplots.pdf",
        path_plot_hyperparams="results/ML/{dataset_id}/per_chrom/{chrom}/plots/hyperparams_table.pdf"
    input:
        check="results/ML/{dataset_id}/checks/{chrom}_train_models"
    params:
        performances_dir="results/ML/{dataset_id}/per_chrom/{chrom}/performances",
        hyperparams_dir="results/ML/{dataset_id}/per_chrom/{chrom}/hyperparams",
        selfeats_dir="results/ML/{dataset_id}/per_chrom/{chrom}/selfeats"
    log:
        "logs/select_best_model/{dataset_id}_{chrom}.log"
    conda:
        "../../envs/r_roadmap.yml"
    shell:
        """
        Rscript workflow/scripts/select_best_model.R \
        --input_dir_performances={params.performances_dir} \
        --input_dir_hyperparams={params.hyperparams_dir} \
        --input_dir_selfeats={params.selfeats_dir} \
        --output_path_best_model_name={output.best_model_name} \
        --output_path_consensus_hyperparams={output.consensus_hyperparams} \
        --output_path_consensus_features={output.consensus_features} \
        --output_path_plot_selfeats={output.path_plot_selfeats} \
        --output_path_plot_performance={output.path_plot_performance} \
        --output_path_median_boxplots={output.path_median_boxplots} \
        --output_path_plot_hyperparams={output.path_plot_hyperparams} 2>{log}
        """