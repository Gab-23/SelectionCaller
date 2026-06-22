rule process_bin_level_essentiality:
    output:
        "data/num_overlaps/essential/{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}_{assembly_code}_{window_mode}_{window_size}_{step}"
    input:
        "resources/DepMap/scores/{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}_{assembly_code}_{window_mode}_{window_size}_{step}"
    log:
        "logs/process_bin_level_essentiality/{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}_{assembly_code}_{window_mode}_{window_size}_{step}.log"
    conda:
        "../../envs/r_roadmap.yml"
    shell:
        """
        Rscript workflow/scripts/process_bin_level_essentiality.R \
        --input_path={input} \
        --output_path={output} 2> {log}
        """