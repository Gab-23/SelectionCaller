rule filter_classification_COSMIC:
    output:
        path_phenotype_id="resources/COSMIC/classification/{COSMIC_parameters}"
    input:
        config["COSMIC_classification_locations"][0]
    conda:
        "../../envs/r_backbone.yml"
    log:
        "logs/filter_classification_COSMIC/{COSMIC_parameters}.log"
    shell:
        """
        Rscript workflow/scripts/filter_classification_COSMIC.R \
            --input_path={input} \
            --output_path_phenotype_id={output.path_phenotype_id} \
            --COSMIC_parameters={wildcards.COSMIC_parameters} 2>{log}
        """