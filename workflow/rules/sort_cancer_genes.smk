rule sort_cancer_genes:
    output:
        "resources/IntOGen/processed/{assembly_code}_v{gencode_version}_{tss_subset}_{tumor_type_IntOGen}.sorted.bed"
    input:
        "resources/IntOGen/processed/{assembly_code}_v{gencode_version}_{tss_subset}_{tumor_type_IntOGen}.bed"
    conda:
        "../envs/bedtools.yml"
    log:
        "logs/sort_cancer_genes/{assembly_code}_v{gencode_version}_{tss_subset}_{tumor_type_IntOGen}.log"
    shell:
        """
	     bedtools sort -i {input} > {output} 2> {log}
        """