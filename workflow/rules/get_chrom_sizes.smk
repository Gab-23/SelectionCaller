rule get_chrom_sizes:
    output:
	    "resources/{assembly_code}.chrom.sizes"
    shell:
        """
	    wget -qP resources/ https://hgdownload.cse.ucsc.edu/goldenpath/{wildcards.assembly_code}/bigZips/{wildcards.assembly_code}.chrom.sizes 
        """