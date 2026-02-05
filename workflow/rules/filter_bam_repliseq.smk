rule filter_bam_repliseq:
    output:
        "resources/repliseq/files/bam/{SRA_id}/{assembly_code}.canon.bam"
    input:
        bam="resources/repliseq/files/bam/{SRA_id}/{assembly_code}.sorted.bam",
        index="resources/repliseq/files/bam/{SRA_id}/{assembly_code}.sorted.bam.bai"
    log:
        "logs/filter_bam_repliseq/{assembly_code}_{SRA_id}.log"
    params:
        num_cores_repliseq=num_cores_repliseq
    conda:
        "../envs/bwa_samtools.yml"
    shell:
        """
        samtools view -@ {params.num_cores_repliseq} -o {output} {input.bam} $(echo chr{{1..22}} chrX chrY) 2>{log}
        """

