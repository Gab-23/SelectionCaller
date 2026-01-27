rule filter_CRISPR_models:
    output:
        "resources/DepMap/models/{oncotree_code_DepMap}"
    input:
        models=config["DepMap_models_locations"][0]
    conda:
        "../envs/r_roadmap.yml"
    log:
        "logs/filter_CRISPR_models/{oncotree_code_DepMap}"
    shell:
        """
        Rscript workflow/scripts/filter_CRISPR_models.R \
        --input_path={input} \
        --output_path={output} \
        --oncotree_code_DepMap={wildcards.oncotree_code_DepMap} 2>{log}
        """