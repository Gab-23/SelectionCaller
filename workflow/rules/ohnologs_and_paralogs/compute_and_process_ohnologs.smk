rule compute_and_process_ohnologs:
    output:
        "data/ohnologs_and_paralogs/ohnologs/{assembly_code}_{window_mode}_{window_size}_{step}"
    input:
        path_backbone="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.mask.filtered.sorted.annotated.bed",
        path_ohnologs=ohnologs_location
    log:
        "logs/compute_ohnologs/{assembly_code}_{window_mode}_{window_size}_{step}.log"
    conda:
        "../envs/r_limma.yml"
    shell:
        """
        Rscript workflow/scripts/compute_and_process_ohnologs.R \
        --input_path_backbone={input.path_backbone} \
        --input_path_ohnologs={input.path_ohnologs} \
        --output_path_ohnologs={output}
        """