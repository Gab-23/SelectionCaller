rule filter_backbone:
    output:
	    "resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.filtered.tsv"
    input:
        backbone="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.tsv",
        blacklist="resources/blacklist/{assembly_code}.bed",
        chrom_sizes="resources/reference_data/chrom_sizes/{assembly_code}.clean.chrom.sizes"
    log:
        "logs/filter_backbone/{assembly_code}_{window_mode}_{window_size}_{step}.log"
    conda:
        "../envs/filter_backbone.yml"
    shell:
        """
        bedtools intersect \
        -a {input.backbone} \
        -b {input.blacklist} \
        -g {input.chrom_sizes} \
        -sorted \
        -v > \
        {output} 2>{log}
        """