rule align_reads_hic:
    output:
        "resources/hic/files/bam/{SRA_id}/{assembly_code}.bam"
    input:
        fasta="resources/reference_data/genome/{assembly_code}.fa",
        index="resources/reference_data/genome/{assembly_code}.fa.bwt",
        FWD="resources/hic/files/fastq/{SRA_id}/{SRA_id}_1.fastq.gz",
        REV="resources/hic/files/fastq/{SRA_id}/{SRA_id}_2.fastq.gz"
    log:
        "logs/align_reads_hic/{assembly_code}_{SRA_id}.log"
    params:
        num_cores_hic=num_cores_hic
    conda:
        "../../envs/hic_env.yaml"
    shell:
        """
        bwa mem -t {params.num_cores_hic} -5SP {input.fasta} {input.FWD} {input.REV} | \
        samblaster | \
        samtools view -S -h -b -F 2316 > {output} 2>{log}
        """