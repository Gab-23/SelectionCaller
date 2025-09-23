rule filter_backbone:
    output:
	    "resources/{assembly_code}.{window_mode}.{window_size}.{step}.filtered.backbone.tsv"
    input:
        "resources/{assembly_code}.{window_mode}.{window_size}.{step}.backbone.tsv",
        "resources/{assembly_code}.blacklist.bed",
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