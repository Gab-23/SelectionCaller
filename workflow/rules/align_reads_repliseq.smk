rule align_reads_repliseq:
    output:
        temp("resources/repliseq/files/bam/{SRA_id}/{assembly_code}.bam")
    input:
        FWD="resources/repliseq/files/fastq/{SRA_id}/{SRA_id}_1_val_1.fq.gz",
        REV="resources/repliseq/files/fastq/{SRA_id}/{SRA_id}_2_val_2.fq.gz",
        FASTA="resources/reference_data/genome/{assembly_code}.fa",
        index="resources/reference_data/genome/{assembly_code}.fa.bwt"
    log:
        "logs/align_reads_repliseq/{assembly_code}_{SRA_id}.log"
    params:
        num_cores_repliseq=num_cores_repliseq
    conda:
        "../envs/bwa_samtools.yml"
    shell:
        """
        bwa mem \
        -R '@RG\\tID:{wildcards.SRA_id}\\tSM:{wildcards.SRA_id}\\tLB:lib1\\tPL:ILLUMINA\\tPU:unit1' \
        -t {params.num_cores_repliseq} \
        {input.FASTA} \
        {input.FWD} \
        {input.REV} | samtools view -S -h -b > {output} 2>{log}
        """


