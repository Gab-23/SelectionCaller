rule get_backbone_file:
    output:
	    "resources/{assembly_code}_{window_mode}_{window_size}_{step}_backbone.tsv"
    input:
        "resources/{assembly_code}.chrom.sizes"
    conda:
        "../envs/r_backbone.yml"
    shell:
        """
        Rscript workflow/scripts/get_backbone_file.R \
        --input_path={input} \
        --output_path={output} \
        --window_mode="{wildcards.window_mode}" \
        --window_size={wildcards.window_size} \
        --step={wildcards.step} 
        """