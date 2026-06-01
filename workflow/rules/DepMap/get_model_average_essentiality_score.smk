rule get_model_average_essentiality_score:
    output:
        "resources/DepMap/scores/{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}.bed"
    input:
        path_scores="resources/DepMap/scores/{oncotree_code_DepMap}_{gender_DepMap}",
        path_genes="resources/DepMap/essential/{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}.sorted.bed"
    log:
        "logs/get_model_average_essentiality_score/{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}.log"
    conda:
        "../envs/r_roadmap.yml"
    shell:
        """
        Rscript workflow/scripts/get_model_average_essentiality_score.R \
            --input_path_scores={input.path_scores} \
            --input_path_genes={input.path_genes} \
            --output_path={output} 2> {log}
        """