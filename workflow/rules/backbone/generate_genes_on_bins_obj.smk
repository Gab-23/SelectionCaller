rule generate_genes_on_bins_obj:
    output:
        "resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}.genes_on_bins.RData"
    input:
        "resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}.mask.filtered.sorted.annotated.bed"
    log:
        "logs/generate_genes_on_bins_obj/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}.log"
    conda:
        # NOTE (fix): this rule had no conda: directive at all - under --use-conda, Snakemake
        # falls back to whatever's already on PATH (the launcher's own "snakemake-env"), which
        # doesn't have R installed. Confirmed directly: exit status 127 ("command not found")
        # on a real cluster run. The script only needs R.utils::commandArgs() plus base R
        # (read.table/save) - r_backbone.yml (already used by other rules in this same
        # backbone module) covers this.
        "../../envs/r_backbone.yml"
    shell:
        """
        Rscript workflow/scripts/generate_genes_on_bins_obj.R \
        --input_path={input} \
        --window_size={wildcards.window_size} \
        --output_path={output} 2>{log}
        """
