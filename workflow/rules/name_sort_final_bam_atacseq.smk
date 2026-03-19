rule name_sort_final_bam_atacseq:
    output:
        temp("resources/atacseq/files/bam/{SRA_id}/{assembly_code}_{mapq_thr}.name.sorted.bam")
    input:
        "resources/atacseq/files/bam/{SRA_id}/{assembly_code}_{mapq_thr}.fixed.bam"
    log:
        "logs/name_sort_final_bam_atacseq/{assembly_code}_{mapq_thr}_{SRA_id}.log"
    params:
        num_cores_atacseq=num_cores_atacseq
    conda:
        "../envs/bowtie2_samtools.yml"
    shell:
        """
        samtools sort -@ {params.num_cores_atacseq} -n {input} -o {output} 2>{log}
        """