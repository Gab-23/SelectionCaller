rule filter_coverage_repliseq:
    output:
        path_filtered="resources/repliseq/files/coverage/{SRA_id}/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}.filtered.bed",
        path_plot="results/plots/repliseq/{SRA_id}/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}.pdf"
    input:
        "resources/repliseq/files/coverage/{SRA_id}/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}.bed"
    log:
        "logs/filter_coverage_repliseq/{SRA_id}_{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}.log"
    conda:
        "../../envs/r_plotting.yml"
    shell:
        """
        Rscript workflow/scripts/filter_coverage_repliseq.R \
        --input_path={input} \
        --output_path_plot={output.path_plot} \
        --output_path_filtered={output.path_filtered} \
        --cov_thr={wildcards.cov_thr} 2>{log}
        """