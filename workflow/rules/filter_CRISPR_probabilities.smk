rule filter_CRISPR_probabilities:
    output:
        "resources/DepMap/probabilities/{oncotree_code_DepMap}"
    input:
        IDs = "resources/DepMap/models/{oncotree_code_DepMap}",
        probabilities = config["DepMap_probabilities_locations"][0]
    log:
        "logs/filter_CRISPR_probabilities/{oncotree_code_DepMap}"
    shell:
        """
        head -n 1 {input.probabilities} > {output}
        grep -Fw -f {input.IDs} {input.probabilities} >> {output} 2> {log}
        """