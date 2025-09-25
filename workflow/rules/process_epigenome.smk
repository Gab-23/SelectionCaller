rule process_epigenome:
    output:
	    "data/{eID}_25_imputed12marks{liftover_str}_{file_format}_{assembly_code}_{window_mode}_{window_size}_{step}.tsv"
    input:
        epigenome_path=get_epigenome_bed,
        backbone_path="resources/{assembly_code}_{window_mode}_{window_size}_{step}_filtered_backbone.tsv"
    conda:
        "../envs/r_roadmap.yml"
    shell:
        """
        Rscript workflow/scripts/compute_state_length.R \
        --epigenome_path={input.epigenome_path} \
        --backbone_path={input.backbone_path} \
        --output_path={output}
        """