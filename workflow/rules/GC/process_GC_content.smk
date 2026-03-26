rule process_GC_content:
    output:
        "data/GC/{assembly_code}_{window_mode}_{window_size}_{step}"
    input:
        "resources/GC/{assembly_code}_{window_mode}_{window_size}_{step}"
    conda:
        "../envs/r_backbone.yml"
    log:
        "logs/process_GC_content/{assembly_code}_{window_mode}_{window_size}_{step}"
    shell:
        """
        Rscript workflow/scripts/process_GC_content.R \
        --input_path={input} \
        --output_path={output} 2>{log}
        """