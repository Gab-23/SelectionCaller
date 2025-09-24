rule filter_backbone:
    output:
	    "resources/{assembly_code}_{window_mode}_{window_size}_{step}_filtered_backbone.tsv"
    input:
        "resources/{assembly_code}_{window_mode}_{window_size}_{step}_backbone.tsv",
        "resources/{assembly_code}_blacklist.bed",
        "resources/{assembly_code}.chrom.sizes"
    conda:
        "../envs/filter_backbone.yml"
    shell:
        """
        bedtools intersect \
        -a {input[0]} \
        -b {input[1]} \
        -g {input[2]} \
        -sorted \
        -v > \
        {output};
        rm {input[0]}
        """