rule process_gene_density:
    output:
        essential="resources/num_overlaps/density/essential/{assembly_code}_{window_mode}_{window_size}_{step}_{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}_{tumor_type_IntOGen}",
        OGs="resources/num_overlaps/density/OG/{assembly_code}_{window_mode}_{window_size}_{step}_{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}_{tumor_type_IntOGen}",
        TSGs="resources/num_overlaps/density/TSG/{assembly_code}_{window_mode}_{window_size}_{step}_{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}_{tumor_type_IntOGen}"
    input:
        essential="resources/num_overlaps/density/essential/{assembly_code}_{window_mode}_{window_size}_{step}_{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}.bed",
        OGs="resources/num_overlaps/density/OG/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}_{tumor_type_IntOGen}.bed",
        TSGs="resources/num_overlaps/density/TSG/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}_{tumor_type_IntOGen}.bed"
    log:
        "logs/process_gene_density/{assembly_code}_{window_mode}_{window_size}_{step}_{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}_{tumor_type_IntOGen}.log"
    conda:
        "../envs/r_spdep.yml"
    shell:
        """
        Rscript process_gene_density.R \
        --input_path_essential={input.essential} \
        --input_path_OGs={input.OGs} \
        --input_path_TSGs={input.TSGs} \
        --output_path_essential={output.essential} \
        --output_path_OGs={output.OGs} \
        --output_path_TSGs={output.TSGs}
        """