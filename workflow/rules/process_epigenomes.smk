rule process_epigenome:
    output:
	    "resources/epigenomes/checks/processing/{tissue_type_Roadmap}_{sample_type}_{liftover_str}_{file_format}_{assembly_code}_{window_mode}_{window_size}_{step}"
    input:
        epigenomes_paths="resources/epigenomes/checks/downloads/{tissue_type_Roadmap}_{sample_type}_{liftover_str}_{file_format}",
        backbone_path="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.filtered.sorted.bed"
    conda:
        "../envs/r_roadmap.yml"
    log:
        "logs/process_epigenome/{tissue_type_Roadmap}_{sample_type}_{liftover_str}_{file_format}_{assembly_code}_{window_mode}_{window_size}_{step}.log"
    shell:
        """
        Rscript workflow/scripts/process_epigenomes.R \
        --epigenomes_paths={input.epigenomes_paths} \
        --backbone_path={input.backbone_path} \
        --output_path={output} \
        --is_test_call=F
        """