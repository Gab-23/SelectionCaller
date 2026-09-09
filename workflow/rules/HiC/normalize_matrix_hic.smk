rule normalize_matrix_hic:
    output:
        matrix="resources/hic/files/cooler/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}.balanced.norm.mat",
        bins="resources/hic/files/cooler/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}.balanced.norm.bins"
    input:
        "resources/hic/files/cooler/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}.balanced.cool"
    log:
        "logs/normalize_matrix_hic/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}.log"
    conda:
        "../../envs/hic_processing_env.yml"
    shell:
        """
        python workflow/scripts/normalize_matrix_hic.py \
        --input_path_cool={input} \
        --output_path_matrix={output.matrix} \
        --output_path_bins={output.bins} 2>{log}
        """