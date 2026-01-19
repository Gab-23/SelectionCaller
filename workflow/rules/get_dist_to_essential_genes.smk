rule get_dist_to_essential_genes:
    output:
        "resources/distances/essential/{assembly_code}_{window_mode}_{window_size}_{step}_{oncotree_code_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}.bed"
    input:
        backbone="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.filtered.sorted.bed",
        essential_genes="resources/DepMap/essential/{oncotree_code_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}.sorted.bed"
    conda:
        "../envs/bedtools.yml"
    log:
        "logs/get_dist_to_essential_genes/{assembly_code}_{window_mode}_{window_size}_{step}_{oncotree_code_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}.log"
    shell:
        """
        bedtools closest -a {input.backbone} -b {input.essential_genes} -d | cut -f 1,2,3,4,9 > {output} 2> {log}
        """