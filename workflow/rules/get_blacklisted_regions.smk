rule get_blacklisted_regions:
    output:
	    "resources/blacklist/{assembly_code}.bed"
    input:
        "resources/reference_data/chrom_sizes/{assembly_code}.chrom.sizes"
    log:
        "logs/get_blacklisted_regions/{assembly_code}.log"
    conda:
        "../envs/filter_backbone.yml"
    shell:
        """
	    wget -qO- https://github.com/Boyle-Lab/Blacklist/raw/refs/heads/master/lists/{wildcards.assembly_code}-blacklist.v2.bed.gz | \
        gunzip | \
        bedtools sort -g {input} > {output}
        """