rule get_essential_genes:
    output:
        "resources/DepMap/essential/{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}"
    input:
        probabilities = "resources/DepMap/probabilities/{oncotree_code_DepMap}_{gender_DepMap}",
        scores        = config["DepMap_scores_locations"][0]
    conda:
        "../../envs/r_roadmap.yml"
    log:
        "logs/get_essential_genes/{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}.log"
    shell:
        """
        Rscript workflow/scripts/get_essential_genes.R \
        --input_path={input.probabilities} \
        --scores_path={input.scores} \
        --output_path={output} \
        --essentiality_cutoff={wildcards.essentiality_cutoff} \
        --num_assays_cutoff={wildcards.num_assays_cutoff} 2>{log}
        """