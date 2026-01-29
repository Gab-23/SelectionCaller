rule find_coord_essential_genes:
    output:
        "resources/DepMap/essential/{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}.bed"
    input:
        essential_IDs = "resources/DepMap/essential/{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}",
        protein_coding_genes = "resources/reference_data/gencode/PC/v{gencode_version}_gene.bed"
    log:
        "logs/find_coord_essential_genes/{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}"
    shell:
        """
        grep -Fw -f {input.essential_IDs} {input.protein_coding_genes} | cut -f 1,2,3,5 > {output} 2>{log}
        """