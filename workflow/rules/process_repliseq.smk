rule process_repliseq:
    output:
        path_data="data/repliseq/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_{smoothing_size}",
        path_plot_distribution="results/plots/repliseq/distribution/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_{smoothing_size}.pdf",
        path_plot_landscape="results/plots/repliseq/landscape/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_{smoothing_size}.pdf",
        path_plot_raw_smooth="results/plots/repliseq/raw_smoothed/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_{smoothing_size}.pdf"
    input:
        path_dir="resources/repliseq/files/coverage/",
        path_samples=config["repliseq_selected_samples_locations"][0],
        check="resources/repliseq/checks/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_coverage_filtering"
    log:
        "logs/process_repliseq/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_{smoothing_size}.log"
    conda:
        "../envs/r_plotting.yml"
    shell:
        """
        Rscript workflow/scripts/process_repliseq.R \
        --input_path_dir={input.path_dir} \
        --input_path_samples={input.path_samples} \
        --output_path_data={output.path_data} \
        --output_path_plot_distribution={output.path_plot_distribution} \
        --output_path_plot_landscape={output.path_plot_landscape} \
        --output_path_plot_raw_smooth={output.path_plot_raw_smooth} \
        --window_size={wildcards.window_size} \
        --smoothing_size={wildcards.smoothing_size}
        """