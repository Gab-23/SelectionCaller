rule index_noblacklist_bam_atacseq:
    output:
        temp("resources/atacseq/files/bam/{SRA_id}/{assembly_code}_{mapq_thr}.noblack.sorted.bam.bai")
    input:
        "resources/atacseq/files/bam/{SRA_id}/{assembly_code}_{mapq_thr}.noblack.sorted.bam"
    log:
        "logs/index_bam_atacseq/{assembly_code}_{mapq_thr}_{SRA_id}.log"
    params:
        num_cores_atacseq=num_cores_atacseq
    conda:
        "../envs/bwa_samtools.yml"
    shell:
        """
        samtools index -@ {params.num_cores_atacseq} {input} 2>{log}
        """
