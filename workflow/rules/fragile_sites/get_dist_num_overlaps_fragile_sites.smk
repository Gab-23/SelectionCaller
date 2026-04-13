rule get_dist_num_overlaps_fragile_sites:
    output:
	dist="resources/distances/fragile_sites/{assembly_code}_{window_mode}_{window_size}_{step}.bed",
	num_overlaps="resources/num_overlaps/fragile_sites/{assembly_code}_{window_mode}_{window_size}_{step}.bed"
    input:
        path_fragile_sites="resources/fragile_sites/fragile_sites_li_{assembly_code}.bed",
        path_backbone="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.filtered.sorted.bed"
    log:
        "logs/get_dist_num_overlaps_fragile_sites/{assembly_code}_{window_mode}_{window_size}_{step}.log"
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        bedtools intersect -a {input.path_backbone} -b {input.path_fragile_sites} -c > {output.num_overlaps} 2>> {log}
        bedtools closest -a {input.path_backbone} -b {input.path_fragile_sites} -t first -d | cut -f 1,2,3,4,8 > {output.dist} 2>> {log}
        """
