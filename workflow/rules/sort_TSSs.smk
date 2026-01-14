rule sort_TSSs:
    output:
        "resources/reference_data/gencode/{assembly_code}_v{gencode_version}_{tss_subset}.sorted.bed"
    input:
        "resources/reference_data/gencode/{assembly_code}_v{gencode_version}_{tss_subset}.bed"
    conda:
        "../envs/bedtools.yml"
    log:
        "logs/sort_TSSs/{assembly_code}_v{gencode_version}_{tss_subset}.log"
    shell:
        """
	     bedtools sort -i {input} > {output} 2> {log}
        """