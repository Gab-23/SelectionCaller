rule get_GDC_CNA_data:
    output:
        "resources/GDC/CNAs/checks/{tumor_type_GDC}_{gender_GDC}"
    input:
        "resources/GDC/CNAs/{tumor_type_GDC}_{gender_GDC}_manifest"
    log:
        "logs/get_GDC_CNA_data/{tumor_type_GDC}_{gender_GDC}.log"
    conda:
        "../envs/gdc.yml"
    shell:
        """
        touch {output}
        gdc-client download -m {input} -d resources/GDC/CNAs/files/
        """