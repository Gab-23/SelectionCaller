rule filter_deduplicated_bam_atacseq:
    output:
        temp("resources/atacseq/files/bam/{SRA_id}/{assembly_code}_{mapq_thr}.dedup.filtered.bam")
    input:
        "resources/atacseq/files/bam/{SRA_id}/{assembly_code}.dedup.bam",
    log:
        "logs/filter_bam_atacseq/{assembly_code}_{mapq_thr}_{SRA_id}.log"
    params:
        num_cores_atacseq=num_cores_atacseq
    conda:
        "../envs/bowtie2_samtools.yml"
    shell:
        """
        samtools view -@ {params.num_cores_atacseq} \
        -q {wildcards.mapq_thr} \
        -f 0x2 \
        -F 0x4 \
        -F 0x8 \
        -F 0x100 \
        -F 0x200 \
        -F 0x400 \
        -F 0x800 \
        -b {input} > {output} 2>{log}
        """

