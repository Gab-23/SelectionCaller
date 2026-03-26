rule get_samples_GDC:
    output:
        path_table="resources/GDC/CNAs/{tumor_type_GDC}_{gender_GDC}.tsv",
        path_barcodes="resources/GDC/CNAs/{tumor_type_GDC}_{gender_GDC}_barcodes",
        path_files="resources/GDC/CNAs/{tumor_type_GDC}_{gender_GDC}_files"
    input:
        path_clinical=config["GDC_clinical_locations"][0],
        path_sample_sheet=config["GDC_sample_sheet_locations"][0]
    params:
        cohort_GDC=cohort_GDC
    conda:
        "../envs/r_backbone.yml"
    log:
        "logs/get_samples_table/{tumor_type_GDC}_{gender_GDC}.log"
    shell:
        """
        Rscript workflow/scripts/get_samples_GDC.R \
        --input_path_clinical={input.path_clinical} \
        --input_path_sample_sheet={input.path_sample_sheet} \
        --gender_GDC={wildcards.gender_GDC} \
        --cohort_GDC={params.cohort_GDC} \
        --output_path_table={output.path_table} \
        --output_path_barcodes={output.path_barcodes} \
        --output_path_files={output.path_files} 2>{log}
        """
    