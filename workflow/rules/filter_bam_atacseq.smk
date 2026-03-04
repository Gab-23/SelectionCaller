rule filter_bam_atacseq:
    output:
        temp("resources/atacseq/files/bam/{SRA_id}/{assembly_code}.canon.bam")
    input:
        bam="resources/atacseq/files/bam/{SRA_id}/{assembly_code}.sorted.bam",
        index="resources/atacseq/files/bam/{SRA_id}/{assembly_code}.sorted.bam.bai"
    log:
        "logs/filter_bam_atacseq/{assembly_code}_{SRA_id}.log"
    params:
        num_cores_atacseq=num_cores_atacseq
    conda:
        "../envs/bowtie2_samtools.yml"
    shell:
        """
        samtools view -@ {params.num_cores_atacseq} -o {output} {input.bam} $(echo chr{{1..22}} chrX chrY) 2>{log}
        """

