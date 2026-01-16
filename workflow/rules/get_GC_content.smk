rule get_GC_content:
    output:
        "resources/GC/{assembly_code}_{window_mode}_{window_size}_{step}"
    input:
        "resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.filtered.sorted.bed"
    conda:
        "../envs/gc_content.yml"
    log:
        "logs/get_gc_content/{assembly_code}_{window_mode}_{window_size}_{step}.log"
    shell:
        """
        Rscript workflow/scripts/get_GC_content.R \
        --input_path={input} \
        --output_path={output} \
        --assembly_code={wildcards.assembly_code} 2>{log}
        """