rule aggregate_epigenomes:
    output:
        "data/epigenomes/{tissue_type_Roadmap}_{sample_type}_{liftover_str}_{file_format}_{assembly_code}_{window_mode}_{window_size}_{step}"
    input:
        "resources/epigenomes/checks/processing/{tissue_type_Roadmap}_{sample_type}_{liftover_str}_{file_format}_{assembly_code}_{window_mode}_{window_size}_{step}"
    conda:
        "../../envs/r_roadmap.yml"
    log:
        "logs/aggregate_epigenomes/{tissue_type_Roadmap}_{sample_type}_{liftover_str}_{file_format}_{assembly_code}_{window_mode}_{window_size}_{step}.log"
    shell:
        """
        Rscript workflow/scripts/aggregate_epigenomes.R \
        --input_path={input} \
        --output_path={output} 2>{log}
        """
