rule sort_TSSs:
    output:
        "resources/reference_data/gencode/v{gencode_version}.sorted.bed"
    input:
        "resources/reference_data/gencode/v{gencode_version}.bed"
    conda:
        "../envs/bedtools.yml"
    log:
        "logs/sort_TSSs/v{gencode_version}.log"
    shell:
        """
	    bedtools sort -i {input} > {output} 2> {log}
        """