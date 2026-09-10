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
        # NOTE (fix): the glob was "resources/hic/files/cooler/*/*.cool", which grabs every
        # per-sample .cool file under ANY subdirectory, regardless of assembly_code/mapq_thr/
        # window_mode/window_size/step - if this pipeline is ever run with more than one
        # combination of these parameters (each has its own config list, like every other
        # parameter in this codebase), this would merge .cool files from DIFFERENT parameter
        # combinations together. Harmless with the current single-value config, but wrong in
        # general. Fixed to only match this specific combination's files.
        """
        COOLER_PATHS=$(ls resources/hic/files/cooler/*/{wildcards.assembly_code}_{wildcards.mapq_thr}_{wildcards.window_mode}_{wildcards.window_size}_{wildcards.step}.cool) 2>>{log}
        cooler merge {output} $COOLER_PATHS 2>>{log}
        """
