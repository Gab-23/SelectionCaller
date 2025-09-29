rule process_epigenome:
    output:
	    "resources/epigenomes_processing_log_{tissue_type}_{liftover_str}_{file_format}_{assembly_code}_{window_mode}_{window_size}_{step}"
    input:
        epigenomes_paths="resources/roadmap_epigenomes_download_log_{tissue_type}_{liftover_str}_{file_format}",
        backbone_path="resources/{assembly_code}_{window_mode}_{window_size}_{step}_filtered_backbone.tsv"
    conda:
        "../envs/r_roadmap.yml"
    shell:
        """
        Rscript workflow/scripts/process_epigenomes.R \
        --epigenomes_paths={input.epigenomes_paths} \
        --backbone_path={input.backbone_path} \
        --output_path={output}
        """