rule balance_cooler:
    output:
        "resources/hic/files/cooler/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}_{hic_cell_line}.balanced.cool"
    input:
        "resources/hic/files/cooler/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}_{hic_cell_line}.cool"
    log:
        "logs/balance_cooler/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}_{hic_cell_line}.log"
    params:
        num_cores_hic=num_cores_hic
    threads:
        num_cores_hic
    resources:
        mem_mb=24000
    conda:
        "../../envs/hic_processing_env.yml"
    shell:
        """
        cp {input} {output}
        cooler balance -p {params.num_cores_hic} {output} 2>{log}
        """