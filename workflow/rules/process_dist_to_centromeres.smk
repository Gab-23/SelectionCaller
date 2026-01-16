rule process_dist_to_centromeres:
    output:
        "data/distances/centromeres/{assembly_code}_{window_mode}_{window_size}_{step}"
    input:
        "resources/distances/centromeres/{assembly_code}_{window_mode}_{window_size}_{step}.bed"
    log:
        "logs/process_dist_to_centromeres/{assembly_code}_{window_mode}_{window_size}_{step}.log"
    shell:
        """
        Rscript workflow/scripts/process_dist_to_centromeres.R \
            --input_path={input} \
            --output_path={output} 2> {log}
        """