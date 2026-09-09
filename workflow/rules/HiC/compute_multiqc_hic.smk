rule compute_multiqc_hic:
    output:
        directory("resources/hic/files/multiqc/{SRA_id}/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}")
    input:
        "resources/hic/files/pairs/{SRA_id}/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}.filtered.dedup.clipped.stats"
    log:
        "logs/compute_multiqc_hic/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{SRA_id}.log"
    conda:
        "../../envs/hic_processing_env.yml"
    shell:
        """
        multiqc --outdir {output} {input} 2>{log}
        """
