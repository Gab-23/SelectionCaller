rule generate_genes_on_bins_obj:
    output:
        "resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}.genes_on_bins.RData"
    input:
        "resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}.mask.filtered.sorted.annotated.bed"
    log:
        "logs/generate_genes_on_bins_obj/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}.log"
    shell:
        """
        Rscript workflow/scripts/generate_genes_on_bins_obj.R \
        --input_path={input} \
        --window_size={wildcards.window_size} \
        --output_path={output} 2> {log}
        """
