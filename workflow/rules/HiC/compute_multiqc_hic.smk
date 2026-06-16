rule compute_multiqc_hic:
    output:
        directory("resources/hic/files/multiqc/{SRA_id}/{assembly_code}_{mapq_thr}")
    input:
        "resources/hic/files/pairs/{SRA_id}/{assembly_code}_{mapq_thr}.filtered.dedup.stats"
    log:
        "logs/compute_multiqc_hic/{assembly_code}_{mapq_thr}_{SRA_id}.log"
    conda:
        "../envs/hic_processing_env.yaml"
    shell:
        """
        multiqc --outdir {output} {input} 2>{log}
        """
