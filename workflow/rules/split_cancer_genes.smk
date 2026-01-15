rule split_cancer_genes:
    output:
        path_TSG="resources/IntOGen/TSG/{assembly_code}_v{gencode_version}_{tumor_type_IntOGen}.sorted.bed",
        path_OG="resources/IntOGen/OG/{assembly_code}_v{gencode_version}_{tumor_type_IntOGen}.sorted.bed"
    input:
        "resources/IntOGen/processed/{assembly_code}_v{gencode_version}_{tumor_type_IntOGen}.sorted.bed"
    log:
        "logs/split_cancer_genes/{assembly_code}_v{gencode_version}_{tumor_type_IntOGen}.log"
    shell:
        """
	     grep -Fw TSG {input} > {output.path_TSG}
         grep -Fw OG {input} > {output.path_OG}
        """