rule get_dist_closest_TSG:
    output:
        "resources/distances/TSG/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}_{tumor_type_IntOGen}.bed"
    input:
        path_backbone="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.filtered.sorted.bed",
        path_TSG="resources/IntOGen/TSG/{assembly_code}_v{gencode_version}_{tumor_type_IntOGen}.sorted.bed"
    log:
        "logs/get_dist_closest_TSG/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}_{tumor_type_IntOGen}.log"
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        bedtools closest -a {input.path_backbone} -b {input.path_TSG} -d | cut -f 1,2,3,4,9,10,11 > {output} 2>{log}
        """