rule sort_model_average_essentiality_score:
    output:
        "resources/DepMap/scores/{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}.sorted.bed"
    input:
        "resources/DepMap/scores/{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}.bed"
    conda:
        "../../envs/bedtools.yml"
    log:
        "logs/sort_model_average_essentiality_score/{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}.log"
    shell:
        """
	    bedtools sort -i {input} > {output} 2> {log}
        """