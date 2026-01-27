rule process_dist_to_essential_genes:
    output:
        "data/distances/essential/{assembly_code}_{window_mode}_{window_size}_{step}_{oncotree_code_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}.bed"
    input:
        "resources/distances/essential/{assembly_code}_{window_mode}_{window_size}_{step}_{oncotree_code_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}.bed"
    log:
        "logs/process_dist_to_essential_genes/{assembly_code}_{window_mode}_{window_size}_{step}_{oncotree_code_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}.log"
    conda:
        "../envs/r_roadmap.yml"
    shell:
        """
        Rscript workflow/scripts/process_dist_to_essential_genes.R \
            --input_path={input} \
            --output_path={output} 2>{log}
        """