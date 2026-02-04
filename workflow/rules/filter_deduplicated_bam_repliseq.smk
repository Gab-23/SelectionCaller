rule filter_deduplicated_bam_repliseq:
    output:
        "resources/repliseq/files/bam/{assembly_code}_{mapq_thr}_{SRA_id}.dedup.filtered.bam"
    input:
        "resources/repliseq/files/bam/{assembly_code}_{SRA_id}.dedup.bam",
    log:
        "logs/filter_bam_repliseq/{assembly_code}_{mapq_thr}_{SRA_id}.log"
    params:
        num_cores_repliseq=num_cores_repliseq
    conda:
        "../envs/bwa_samtools.yml"
    shell:
        """
        samtools view -@ {params.num_cores_repliseq} -q {wildcards.mapq_thr} -b {input} > {output} 2>{log}
        """

