rule process_repliseq:
    # NOTE (30/10/2026): wildcard {repliseq_cell_line} (linea cellulare del Repli-seq, dalla mappatura tessuto -> linea), nel
    # nome di tutte le uscite: linee diverse (MCF10A per il seno, T2P per il LUAD) non condividono piu' lo stesso file.
    # Il file dei campioni e' quello della linea: in config, selected_samples_location_repliseqs deve contenere il segnaposto
    # {repliseq_cell_line} (es. "config/{repliseq_cell_line}_selected_samples_repliseq").
    # I replicati (colonna facoltativa "replicate" del file dei campioni) vengono combinati sommando le coperture per
    # sample_type prima della CPM; path_plot_replicates contiene i grafici per replicato e combinati (vedi process_repliseq.R).
    # Anche le dimensioni di libreria e il controllo delle coperture filtrate sono specifici per la linea cellulare.
    output:
        path_data="data/repliseq/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_{smoothing_size}_{repliseq_cell_line}",
        path_plot_distribution="results/plots/repliseq/distribution/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_{smoothing_size}_{repliseq_cell_line}.pdf",
        path_plot_landscape="results/plots/repliseq/landscape/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_{smoothing_size}_{repliseq_cell_line}.pdf",
        path_plot_raw_smooth="results/plots/repliseq/raw_smoothed/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_{smoothing_size}_{repliseq_cell_line}.pdf",
        path_plot_replicates="results/plots/repliseq/replicates/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_{smoothing_size}_{repliseq_cell_line}.pdf"
    input:
        path_dir="resources/repliseq/files/coverage/",
        path_samples=selected_samples_location_repliseq,
        path_library_size="resources/repliseq/files/library_size/{assembly_code}_{mapq_thr}_{repliseq_cell_line}",
        check="resources/repliseq/checks/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_{repliseq_cell_line}_coverage_filtering",
        backbone_mask = "resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.mask.filtered.sorted.bed"
    log:
        "logs/process_repliseq/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_{smoothing_size}_{repliseq_cell_line}.log"
    conda:
        "../../envs/r_process_repliseq.yml"
    shell:
        """
        Rscript workflow/scripts/process_repliseq.R \
        --input_path_dir={input.path_dir} \
        --input_path_samples={input.path_samples} \
        --input_path_library_size={input.path_library_size} \
        --input_path_backbone_mask={input.backbone_mask} \
        --output_path_data={output.path_data} \
        --output_path_plot_distribution={output.path_plot_distribution} \
        --output_path_plot_landscape={output.path_plot_landscape} \
        --output_path_plot_raw_smooth={output.path_plot_raw_smooth} \
        --output_path_plot_replicates={output.path_plot_replicates} \
        --window_size={wildcards.window_size} \
        --smoothing_size={wildcards.smoothing_size} \
        --assembly_code={wildcards.assembly_code} \
        --window_mode={wildcards.window_mode} \
        --step={wildcards.step} \
        --mapq_thr={wildcards.mapq_thr} \
        --cov_thr={wildcards.cov_thr} 2>{log}
        """