rule process_dist_to_telomeres:
    output:
        "data/distances/telomeres/{assembly_code}_{window_mode}_{window_size}_{step}"
    input:
        "resources/distances/telomeres/{assembly_code}_{window_mode}_{window_size}_{step}.bed"
    conda:
        "../../envs/r_roadmap.yml"
    log:
        "logs/process_dist_to_telomeres/{assembly_code}_{window_mode}_{window_size}_{step}.log"
    shell:
        """
        Rscript workflow/scripts/process_dist_to_telomeres.R \
            --input_path={input} \
            --output_path={output} 2> {log}
        """