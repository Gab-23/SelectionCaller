rule get_coding_mutations_description:
    output:
        "resources/COSMIC/mutations/coding/{gender_COSMIC}_{COSMIC_parameters}_w_description.bed"
    input:
        path_coords="resources/COSMIC/coordinates/coding/{gender_COSMIC}_{COSMIC_parameters}.sorted.bed",
        path_muts="resources/COSMIC/mutations/coding/{gender_COSMIC}_{COSMIC_parameters}.filtered.bed"
    conda:
        "../envs/r_roadmap.yml"
    log:
        "logs/get_coding_mutations_description/{gender_COSMIC}_{COSMIC_parameters}.log"
    shell:
        """
        Rscript workflow/scripts/get_coding_mutations_description.R \
        --input_path_coords={input.path_coords} \
        --input_path_muts={input.path_muts} \
        --output_path={output} 2>{log}
        """