rule get_essential_genes:
    output:
        "resources/DepMap/essential/{oncotree_code_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}"
    input:
        "resources/DepMap/probabilities/{oncotree_code_DepMap}"
    conda:
        "../envs/r_roadmap.yml"
    log:
        "logs/filter_CRISPR_probabilities/{oncotree_code_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}.log"
    shell:
        """
        Rscript workflow/scripts/get_essential_genes.R \
        --input_path={input} \
        --output_path={output} \
        --essentiality_cutoff={wildcards.essentiality_cutoff} \
        --num_assays_cutoff={wildcards.num_assays_cutoff}
        """