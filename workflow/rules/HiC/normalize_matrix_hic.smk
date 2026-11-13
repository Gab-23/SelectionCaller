rule normalize_matrix_hic:
    output:
        matrix="resources/hic/files/cooler/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}_{hic_cell_line}_{chrom}.balanced.norm.mat",
        bins="resources/hic/files/cooler/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}_{hic_cell_line}_{chrom}.balanced.norm.bins"
    input:
        "resources/hic/files/cooler/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}_{hic_cell_line}.balanced.cool"
    log:
        "logs/normalize_matrix_hic/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}_{hic_cell_line}_{chrom}.log"
    conda:
        "../../envs/hic_processing_env.yml"
    shell:
        """
        python workflow/scripts/normalize_matrix_hic.py \
        --input_path_cool={input} \
        --chrom={wildcards.chrom} \
        --output_path_matrix={output.matrix} \
        --output_path_bins={output.bins} 2>{log}
        """