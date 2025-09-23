rule get_blacklisted_regions:
    output:
	    "resources/{assembly_code}.blacklist.bed"
    input:
        "resources/{assembly_code}.chrom.sizes"
    conda:
        "../envs/filter_backbone.yml"
    shell:
        """
	    wget -qO- https://github.com/Boyle-Lab/Blacklist/raw/refs/heads/master/lists/{wildcards.assembly_code}-blacklist.v2.bed.gz | \
        gunzip -c | \
        bedtools sort -g {input} > {output}
        """