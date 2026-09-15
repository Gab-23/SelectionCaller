rule process_bins_hic:
    output:
        "resources/hic/files/cooler/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}.balanced.norm.processed.bins"
    input:
        path_backbone_raw="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.sorted.bed",
        path_backbone_filtered="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.filtered.sorted.bed",
        path_bins="resources/hic/files/cooler/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}.balanced.norm.bins"
    log:
        "logs/process_bins_hic/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}.log"
    conda:
        "../../envs/r_backbone.yml"
    shell:
        """
        Rscript workflow/scripts/process_bins_hic.R \
        --input_path_backbone_raw={input.path_backbone_raw} \
        --input_path_backbone_filtered={input.path_backbone_filtered} \
        --input_path_bins={input.path_bins} \
        --output_path={output} 2>{log}
        """