rule compute_coding_mutations_overlap:
    output:
        "data/num_overlaps/mutations/coding/{assembly_code}_{window_mode}_{window_size}_{step}_{gender_COSMIC}_{COSMIC_parameters}"
    input:
        path_mutation="resources/COSMIC/mutations/coding/{gender_COSMIC}_{COSMIC_parameters}_w_description.sorted.bed",
        path_backbone="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.mask.filtered.sorted.bed"
    conda:
        "../envs/r_roadmap.yml"
    log:
        "logs/compute_coding_mutations_overlap/{assembly_code}_{window_mode}_{window_size}_{step}_{gender_COSMIC}_{COSMIC_parameters}.log"
    shell:
        """
        Rscript workflow/scripts/compute_coding_mutations_overlap.R \
        --input_path_mutation={input.path_mutation} \
        --input_path_backbone={input.path_backbone} \
        --output_path_mutation={output} 2>{log}
        """