rule process_GTEx:
    # Gene x TCGA-cohort matrix of tissue-relative normal expression (z-scores of median
    # TPM across all mapped cohorts). Gene-level: independent of assembly_code.
    output:
        "data/GTEx/{GTEx_release}/tissue_relative_expression.tsv"
    input:
        gct="resources/GTEx/{GTEx_release}/gene_median_tpm.gct.gz",
        mapping=mapping_file_GTEx_location
    wildcard_constraints:
        GTEx_release="v[0-9]+"
    log:
        "logs/process_GTEx/{GTEx_release}.log"
    conda:
        "../../envs/r_plain.yml"
    shell:
        """
        Rscript workflow/scripts/process_GTEx.R \
        --input_path_gct={input.gct} \
        --input_path_mapping={input.mapping} \
        --output_path={output} 2> {log}
        """