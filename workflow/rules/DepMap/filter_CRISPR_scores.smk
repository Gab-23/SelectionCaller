rule filter_CRISPR_scores:
    output:
        "resources/DepMap/scores/{oncotree_code_DepMap}_{gender_DepMap}"
    input:
        IDs = "resources/DepMap/models/{oncotree_code_DepMap}_{gender_DepMap}",
        scores = config["DepMap_scores_locations"][0]
    log:
        "logs/filter_CRISPR_scores/{oncotree_code_DepMap}_{gender_DepMap}"
    shell:
        """
        head -n 1 {input.scores} > {output}
        grep -Fw -f {input.IDs} {input.scores} >> {output} 2> {log}
        """