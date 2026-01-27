rule filter_classification_COSMIC:
    output:
        path_filtered="resources/COSMIC/classification/filtered_classification.tsv",
        path_phenotype_id="resources/COSMIC/classification/filtered_phenotype_ids.tsv"
    input:
        config["COSMIC_classification_locations"][0]
    params:
        COSMIC_parameters=COSMIC_parameters
    conda:
        "../envs/r_backbone.yml"
    log:
        "logs/filter_classification_COSMIC/filter_classification_COSMIC.log"
    shell:
        """
        Rscript workflow/scripts/filter_classification_COSMIC.R \
            --input_path={input} \
            --output_path_filtered={output.path_filtered} \
            --output_path_phenotype_id={output.path_phenotype_id} \
            --COSMIC_parameters={params.COSMIC_parameters} 2>{log}
        """