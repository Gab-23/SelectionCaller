rule aggregate_epigenomes:
    output:
        "data/{tissue_type}_aggregated_{liftover_str}_{file_format}_{assembly_code}_{window_mode}_{window_size}_{step}.tsv"
    input:
        "resources/epigenomes_processing_log_{tissue_type}_{liftover_str}_{file_format}_{assembly_code}_{window_mode}_{window_size}_{step}"
    conda:
        "../envs/r_roadmap.yml"
    shell:
        """
        Rscript workflow/scripts/aggregate_epigenomes.R \
        --input_path={input} \
        --output_path={output}
        """
