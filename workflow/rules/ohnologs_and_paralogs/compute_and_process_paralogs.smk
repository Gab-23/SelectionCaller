rule compute_and_process_paralogs:
    output:
        "data/ohnologs_and_paralogs/paralogs/{assembly_code}_{window_mode}_{window_size}_{step}"
    input:
        path_backbone="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.mask.filtered.sorted.annotated.bed",
        path_paralogs=paralogs_location
    log:
        "logs/compute_paralogs/{assembly_code}_{window_mode}_{window_size}_{step}.log"
    conda:
        "../../envs/r_limma.yml"
    shell:
        """
        Rscript workflow/scripts/compute_and_process_paralogs.R \
        --input_path_backbone={input.path_backbone} \
        --input_path_paralogs={input.path_paralogs} \
        --output_path_paralogs={output}
        """