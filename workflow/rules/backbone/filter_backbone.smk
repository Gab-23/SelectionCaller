rule filter_backbone:
    output:
	    backbone_full="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.filtered.bed",
        backbone_mask="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.mask.filtered.bed"
    input:
        backbone_full="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.bed",
        backbone_mask="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.mask.bed",
        blacklist="resources/blacklist/{assembly_code}.bed",
        chrom_sizes="resources/reference_data/chrom_sizes/{assembly_code}.chrom.sizes"
    log:
        "logs/filter_backbone/{assembly_code}_{window_mode}_{window_size}_{step}.log"
    conda:
        "../envs/filter_backbone.yml"
    shell:
        """
        bedtools intersect -a {input.backbone_full} -b {input.blacklist} -g {input.chrom_sizes} -sorted -v > {output.backbone_full} 2>>{log}
        bedtools intersect -a {input.backbone_mask} -b {input.blacklist} -g {input.chrom_sizes} -sorted -v > {output.backbone_mask} 2>>{log}
        """