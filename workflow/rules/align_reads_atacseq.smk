rule align_reads_atacseq:
    output:
        "resources/atacseq/files/bam/{SRA_id}/{assembly_code}.bam"
    input:
        FWD="resources/atacseq/files/fastq/{SRA_id}/{SRA_id}_1_val_1.fq.gz",
        REV="resources/atacseq/files/fastq/{SRA_id}/{SRA_id}_2_val_2.fq.gz",
        FASTA="resources/reference_data/genome/{assembly_code}.fa",
        index="resources/reference_data/genome/{assembly_code}.1.bt2"
    log:
        "logs/align_reads_atacseq/{assembly_code}_{SRA_id}.log"
    params:
        num_cores_atacseq=num_cores_atacseq,
        max_insert_size_atacseq=max_insert_size_atacseq,
        index_prefix="resources/reference_data/genome/{assembly_code}"
    conda:
        "../envs/bowtie2_samtools.yml"
    shell:
        """
        bowtie2 -q \
        -p {params.num_cores_atacseq} \
        -X {params.max_insert_size_atacseq} \
        --local \
        --very-sensitive-local \
        --no-discordant \
        --no-mixed \
        --dovetail \
		-x {input.index} \
        -1 {input.FWD} \
        -2 {input.REV} | \
        samtools view -@ {params.num_cores_atacseq} \
        -b \
        -o {output} 2>{log}
        """


