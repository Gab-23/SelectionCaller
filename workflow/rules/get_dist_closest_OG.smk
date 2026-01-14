rule get_dist_closest_OG:
    output:
        "resources/IntOGen/OG/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}_{tss_subset}_{tumor_type_IntOGen}.bed"
    input:
        path_backbone="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.filtered.sorted.bed",
        path_OG="resources/IntOGen/OG/{assembly_code}_v{gencode_version}_{tss_subset}_{tumor_type_IntOGen}.sorted.bed"
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        bedtools closest -a {input.path_backbone} -b {input.path_OG} -d | cut -f 1,2,3,4,9,10,11 > {output}
        """