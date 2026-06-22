rule filter_CRISPR_models:
    output:
        "resources/DepMap/models/{oncotree_code_DepMap}_{gender_DepMap}"
    input:
        models=config["DepMap_models_locations"][0]
    conda:
        "../../envs/r_roadmap.yml"
    log:
        "logs/filter_CRISPR_models/{oncotree_code_DepMap}_{gender_DepMap}"
    shell:
        """
        Rscript workflow/scripts/filter_CRISPR_models.R \
        --input_path={input} \
        --output_path={output} \
        --gender_DepMap={wildcards.gender_DepMap} \
        --oncotree_code_DepMap={wildcards.oncotree_code_DepMap} 2>{log}
        """