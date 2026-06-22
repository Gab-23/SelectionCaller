rule sort_noblacklist_bam_atacseq:
    output:
        temp("resources/atacseq/files/bam/{SRA_id}/{assembly_code}_{mapq_thr}.noblack.sorted.bam")
    input:
        "resources/atacseq/files/bam/{SRA_id}/{assembly_code}_{mapq_thr}.noblack.bam"
    log:
        "logs/sort_bam_atacseq/{assembly_code}_{mapq_thr}_{SRA_id}.log"
    params:
        num_cores_atacseq=num_cores_atacseq
    conda:
        "../../envs/bwa_samtools.yml"
    shell:
        """
        samtools sort -@ {params.num_cores_atacseq} -o {output} {input} 2>{log}
        """
