rule sort_backbone:
    output:
        backbone_full="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.filtered.sorted.bed",
        backbone_mask="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.mask.filtered.sorted.bed",
    input:
        backbone_full="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.filtered.bed",
        backbone_mask="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.mask.filtered.bed"
    conda:
        "../envs/bedtools.yml"
    log:
        "logs/sort_backbone/{assembly_code}_{window_mode}_{window_size}_{step}.log"
    shell:
        """
         bedtools sort -i {input.backbone_full} > {output.backbone_full} 2>> {log}
         bedtools sort -i {input.backbone_mask} > {output.backbone_mask} 2>> {log}
        """