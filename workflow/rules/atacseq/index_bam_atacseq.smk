rule index_bam_atacseq:
    output:
        temp("resources/atacseq/files/bam/{SRA_id}/{assembly_code}.sorted.bam.bai")
    input:
        "resources/atacseq/files/bam/{SRA_id}/{assembly_code}.sorted.bam"
    log:
        "logs/index_bam_atacseq/{assembly_code}_{SRA_id}.log"
    params:
        num_cores_atacseq=num_cores_atacseq
    conda:
        "../../envs/bowtie2_samtools.yml"
    shell:
        """
        samtools index -@ {params.num_cores_atacseq} {input} 2>{log}
        """
