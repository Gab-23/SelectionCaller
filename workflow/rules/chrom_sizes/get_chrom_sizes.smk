rule get_chrom_sizes:
    output:
	    "resources/reference_data/chrom_sizes/{assembly_code}.chrom.sizes"
    log:
        "logs/get_chrom_sizes/{assembly_code}.log"
    shell:
        """
	    wget -qP resources/reference_data/chrom_sizes/ https://hgdownload.cse.ucsc.edu/goldenpath/{wildcards.assembly_code}/bigZips/{wildcards.assembly_code}.chrom.sizes 2> {log}
        """