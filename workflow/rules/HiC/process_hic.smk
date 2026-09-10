rule process_hic:
    output:
        path_eigen="data/hic/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}.eigen"
    input:
        path_matrix="resources/hic/files/cooler/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}.balanced.norm.mat",
        path_bins="resources/hic/files/cooler/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}.balanced.norm.processed.bins"
    log:
        "logs/process_hic/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}.log"
    conda:
        "../../envs/hic_processing_env.yml"
    shell:
        """
        python workflow/scripts/process_hic.py \
        --input_path_matrix={input.path_matrix} \
        --input_path_bins={input.path_bins} \
        --output_path_eigen={output.path_eigen}
        """
