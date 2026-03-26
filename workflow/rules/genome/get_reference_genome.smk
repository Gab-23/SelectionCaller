rule get_reference_genome:
    output:
        "resources/reference_data/genome/{assembly_code}.fa"
    log:
        "logs/get_reference_genome/{assembly_code}.log"
    shell:
        """
        wget -qO- https://hgdownload.cse.ucsc.edu/goldenPath/{assembly_code}/bigZips/latest/{assembly_code}.fa.gz | gunzip > {output}  2>> {log}
        """