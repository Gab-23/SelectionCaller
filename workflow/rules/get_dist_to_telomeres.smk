rule get_dist_to_telomeres:
    output:
        "resources/distances/telomeres/{assembly_code}_{window_mode}_{window_size}_{step}.bed"
    input:
        backbone="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.filtered.sorted.bed",
        telomeres="resources/reference_data/telomeres/{assembly_code}.bed"
    conda:
        "../envs/bedtools.yml"
    log:
        "logs/get_dist_to_telomeres/{assembly_code}_{window_mode}_{window_size}_{step}.log"
    shell:
        """
        bedtools closest -a {input.backbone} -b {input.telomeres} -d | cut -f 1,2,3,4,8 > {output} 2> {log}
        """