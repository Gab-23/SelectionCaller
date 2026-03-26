rule clean_chrom_sizes:
    output:
        "resources/reference_data/chrom_sizes/{assembly_code}.clean.chrom.sizes"
    input:
        "resources/reference_data/chrom_sizes/{assembly_code}.chrom.sizes"
    log:
        "logs/clean_chrom_sizes/{assembly_code}.log"
    shell:
        """
	    cat {input} | grep -v ["_XMY"] > {output} 2> {log}
        """