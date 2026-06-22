rule sort_essential_genes:
    output:
        "resources/DepMap/essential/{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}.sorted.bed"
    input:
        "resources/DepMap/essential/{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}.bed"
    conda:
        "../../envs/bedtools.yml"
    log:
        "logs/sort_essential_genes/{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}.log"
    shell:
        """
	    bedtools sort -i {input} > {output} 2> {log}
        """