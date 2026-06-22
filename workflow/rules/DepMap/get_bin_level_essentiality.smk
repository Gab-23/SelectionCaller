rule get_bin_level_essentiality:
    output:
        "resources/DepMap/scores/{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}_{assembly_code}_{window_mode}_{window_size}_{step}"
    input:
        path_backbone="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.mask.filtered.sorted.annotated.bed",
        path_essentiality="resources/DepMap/scores/{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}.sorted.bed"
    conda:
        "../../envs/bedtools.yml"
    log:
        "logs/get_bin_level_essentiality/.log"
    shell:
        """
	    bedtools intersect -a {input.path_backbone} -b {input.path_essentiality} -wa -wb | cut -f 4,10 2> {log}
        """