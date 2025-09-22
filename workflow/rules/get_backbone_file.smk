rule get_backbone_file:
    output:
	    "resources/{assembly_code}.{window_mode}.{window_size}.{step}_backbone.tsv"
    input:
        "resources/{assembly_code}.chrom.sizes"
    conda:
        "../envs/r_backbone.yml"
    shell:
        """
        Rscript workflow/scripts/get_backbone_file.R \
        --input_file={input} \
        --output_file={output} \
        --window_size={wildcards.window_size} \
        --window_mode="{wildcards.window_mode}" \
        --step={wildcards.step} 
        """