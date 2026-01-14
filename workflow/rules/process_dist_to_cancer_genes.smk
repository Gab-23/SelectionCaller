rule process_dist_to_cancer_genes:
    output:
        path_TSG="data/distances/TSG/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}_{tss_subset}_{tumor_type_IntOGen}",
        path_OG="data/distances/OG/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}_{tss_subset}_{tumor_type_IntOGen}"
    input:
        path_TSG="resources/IntOGen/TSG/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}_{tss_subset}_{tumor_type_IntOGen}.bed",
        path_OG="resources/IntOGen/OG/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}_{tss_subset}_{tumor_type_IntOGen}.bed"
    conda:
        "../envs/r_roadmap.yml"
    shell:
        """
        Rscript workflow/scripts/process_dist_to_cancer_genes.R \
            --input_path_TSG={input.path_TSG} \
            --input_path_OG={input.path_OG} \
            --output_path_TSG={output.path_TSG} \
            --output_path_OG={output.path_OG}
        """