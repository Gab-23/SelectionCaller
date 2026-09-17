rule process_hic:
    output:
        "data/hic/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}_{chrom}.eigen"
    input:
        matrix="resources/hic/files/cooler/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}_{chrom}.balanced.norm.mat",
        bins="resources/hic/files/cooler/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}_{chrom}.balanced.norm.processed.bins"
    log:
        "logs/process_hic/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}_{chrom}.log"
    conda:
        "../../envs/hic_processing_env.yml"
    shell:
        """
        python workflow/scripts/process_hic.py \
        --input_path_matrix={input.matrix} \
        --input_path_bins={input.bins} \
        --output_path_eigen={output} 2>{log}
        """