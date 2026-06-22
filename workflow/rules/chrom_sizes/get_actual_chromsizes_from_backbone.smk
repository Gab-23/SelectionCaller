rule get_actual_chromsizes_from_backbone:
    output:
        "resources/reference_data/chrom_sizes/{assembly_code}_{window_mode}_{window_size}_{step}.actual.chrom.sizes"
    input:
        "resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.filtered.sorted.bed"
    log:
        "logs/get_actual_chromsizes_from_backbone/{assembly_code}_{window_mode}_{window_size}_{step}.log"
    conda:
        "../../envs/r_backbone.yml"
    shell:
        """
        Rscript workflow/scripts/get_actual_chromsizes_from_backbone.R \
        --input_path={input} \
        --output_path={output}
        """