rule index_noblacklist_bam_repliseq:
    output:
        temp("resources/repliseq/files/bam/{SRA_id}/{assembly_code}_{mapq_thr}.noblack.sorted.bam.bai")
    input:
        "resources/repliseq/files/bam/{SRA_id}/{assembly_code}_{mapq_thr}.noblack.sorted.bam"
    log:
        "logs/index_bam_repliseq/{assembly_code}_{mapq_thr}_{SRA_id}.log"
    params:
        num_cores_repliseq=num_cores_repliseq
    conda:
        "../../envs/bwa_samtools.yml"
    shell:
        """
        samtools index -@ {params.num_cores_repliseq} {input} 2>{log}
        """
