rule index_bam_repliseq:
    output:
        "resources/repliseq/files/bam/{assembly_code}_{SRA_id}.sorted.bam.bai"
    input:
        "resources/repliseq/files/bam/{assembly_code}_{SRA_id}.sorted.bam"
    log:
        "logs/index_bam_repliseq/{assembly_code}_{SRA_id}.log"
    params:
        num_cores_repliseq=num_cores_repliseq
    conda:
        "../envs/bwa_samtools.yml"
    shell:
        """
        samtools index -@ {params.num_cores_repliseq} {input} 2>{log}
        """
