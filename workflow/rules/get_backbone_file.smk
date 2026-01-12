rule get_backbone_file:
    output:
	    "resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.tsv"
    input:
        "resources/reference_data/chrom_sizes/{assembly_code}.clean.chrom.sizes"
    log:
        "logs/get_backbone_file/{assembly_code}_{window_mode}_{window_size}_{step}.log"
    conda:
        "../envs/r_backbone.yml"
    shell:
        """
        Rscript workflow/scripts/get_backbone_file.R \
        --input_path={input} \
        --output_path={output} \
        --window_mode="{wildcards.window_mode}" \
        --window_size={wildcards.window_size} \
        --step={wildcards.step} 2>{log}
        """