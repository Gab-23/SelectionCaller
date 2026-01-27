rule filter_samples_COSMIC:
    output:
        "resources/COSMIC/samples/filtered_samples.tsv"
    input:
        config["COSMIC_samples_locations"][0]
    params:
        sample_name_pattern=config["COSMIC_sample_name_patterns"][0]
    conda:
        "../envs/r_backbone.yml"
    log:
        "logs/filter_samples_COSMIC/filter_samples_COSMIC.log"
    shell:
        """
        Rscript workflow/scripts/filter_samples_COSMIC.R \
            --input_path={input} \
            --output_path={output} \
            --sample_name_pattern={params.sample_name_pattern} 2>{log}
        """