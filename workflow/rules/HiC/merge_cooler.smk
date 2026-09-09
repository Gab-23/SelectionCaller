rule merge_coolers:
    output:
        "resources/hic/files/cooler/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}.cool"
    input:
        check="resources/hic/checks/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}_matrices"
    log:
        "logs/merge_coolers/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}.log"
    conda:
        "../../envs/hic_processing_env.yml"
    shell:
        """
        COOLER_PATHS=$(ls resources/hic/files/cooler/*/*.cool) 2>>{log}
        cooler merge {output} $COOLER_PATHS 2>>{log}
        """