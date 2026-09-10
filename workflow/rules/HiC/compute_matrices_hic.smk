rule compute_matrices_hic:
    output:
        "resources/hic/files/cooler/{SRA_id}/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}.cool"
    input:
        backbone="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.filtered.sorted.bed",
        pairs="resources/hic/files/pairs/{SRA_id}/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}.filtered.dedup.sorted.clipped.pairs"        
    log:
        "logs/compute_matrices_hic/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}_{SRA_id}.log"
    conda:
        "../../envs/hic_processing_env.yml"
    shell:
        """
        cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 \
        --no-symmetric-upper \
        --assembly {wildcards.assembly_code} \
        {input.backbone} \
        {input.pairs} \
        {output} 2>{log}
        """