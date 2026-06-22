rule process_num_overlaps_fragile_sites:
    output:
        "data/num_overlaps/fragile_sites/{assembly_code}_{window_mode}_{window_size}_{step}"
    input:
        "resources/num_overlaps/fragile_sites/{assembly_code}_{window_mode}_{window_size}_{step}.bed",
    log:
        "logs/process_num_overlaps_fragile_sites/{assembly_code}_{window_mode}_{window_size}_{step}.log"
    conda:
        "../../envs/r_roadmap.yml"
    shell:
        """
        Rscript workflow/scripts/process_num_overlaps_fragile_sites.R \
            --input_path_num_overlaps_fragile_sites={input} \
            --output_path_num_overlaps_fragile_sites={output} \
        """
