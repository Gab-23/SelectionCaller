rule get_final_samples_COSMIC:
    output:
        path_samples="resources/COSMIC/samples/{gender_COSMIC}_{COSMIC_parameters}"
    input:
        path_phenotype_ids="resources/COSMIC/classification/{COSMIC_parameters}",
        path_samples="resources/COSMIC/samples/{gender_COSMIC}"
    conda:
        "../../envs/r_backbone.yml"
    log:
        "logs/get_final_samples_COSMIC/{gender_COSMIC}_{COSMIC_parameters}.log"
    shell:
        """
        Rscript workflow/scripts/get_final_samples_COSMIC.R \
            --input_path_samples={input.path_samples} \
            --input_path_phenotype_ids={input.path_phenotype_ids} \
            --output_path_samples={output.path_samples} 2>{log}
        """