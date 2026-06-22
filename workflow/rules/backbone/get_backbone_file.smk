rule get_backbone_file:
    output:
        path_full="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.bed",
        path_mask="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.mask.bed"
    input:
        path_chrom_sizes="resources/reference_data/chrom_sizes/{assembly_code}.clean.chrom.sizes",
        path_mask=""
    log:
        "logs/get_backbone_file/{assembly_code}_{window_mode}_{window_size}_{step}.log"
    conda:
        "../../envs/r_backbone.yml"
    shell:
        """
        Rscript workflow/scripts/get_backbone_file.R \
        --input_path_chrom_sizes={input.path_chrom_sizes} \
        --input_path_mask={input.path_mask} \
        --output_path_full={output.path_full} \
        --output_path_mask={output.path_mask} \
        --window_mode="{wildcards.window_mode}" \
        --window_size={wildcards.window_size} \
        --step={wildcards.step} 2>{log}
        """