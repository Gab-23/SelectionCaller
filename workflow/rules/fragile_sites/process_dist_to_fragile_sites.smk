rule process_dist_to_fragile_sites:
    output:
        "data/distances/fragile_sites/{assembly_code}_{window_mode}_{window_size}_{step}"
    input:
        "resources/distances/fragile_sites/{assembly_code}_{window_mode}_{window_size}_{step}.bed"
    log:
        "logs/process_dist_to_fragile_sites/{assembly_code}_{window_mode}_{window_size}_{step}.log"
    conda:
        "../../envs/r_roadmap.yml"
    shell:
        """
        Rscript workflow/scripts/process_dist_to_fragile_sites.R \
            --input_path_dist_fragile_sites={input} \
            --output_path_dist_fragile_sites={output}
        """
