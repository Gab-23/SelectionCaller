rule sort_filtered_bam_atacseq:
    output:
        temp("resources/atacseq/files/bam/{SRA_id}/{assembly_code}.canon.sorted.bam")
    input:
        "resources/atacseq/files/bam/{SRA_id}/{assembly_code}.canon.bam"
    log:
        "logs/sort_filtered_bam_atacseq/{assembly_code}_{SRA_id}.log"
    params:
        num_cores_atacseq=num_cores_atacseq
    conda:
        "../envs/bowtie2_samtools.yml"
    shell:
        """
        samtools sort -@ {params.num_cores_atacseq} -o {output} {input} 2>{log}
        """
