rule get_dist_to_centromeres:
    output:
        "resources/distances/centromeres/{assembly_code}_{window_mode}_{window_size}_{step}.bed"
    input:
        backbone="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.filtered.sorted.bed",
        centromeres="resources/reference_data/centromeres/{assembly_code}.bed"
    conda:
        "../envs/bedtools.yml"
    log:
        "logs/get_dist_to_centromeres/{assembly_code}_{window_mode}_{window_size}_{step}.log"
    shell:
        """
        bedtools closest -a {input.backbone} -b {input.centromeres} -d | cut -f 1,2,3,4,8 > {output} 2> {log}
        """