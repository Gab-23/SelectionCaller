rule balance_cooler:
    output:
        "resources/hic/files/cooler/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}.balanced.cool"
    input:
        "resources/hic/files/cooler/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}.cool"
    log:
        "logs/balance_cooler/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}.log"
    params:
        num_cores_hic=num_cores_hic
    conda:
        "../../envs/hic_processing_env.yml"
    shell:
        """
        cp {input} {output}
        cooler balance -p {params.num_cores_hic} {output} 2>{log}
        """