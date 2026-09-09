rule trim_reads_atacseq:
    output:
        FWD=temp("resources/atacseq/files/fastq/{SRA_id}/{SRA_id}_1_val_1.fq.gz"),
        REV=temp("resources/atacseq/files/fastq/{SRA_id}/{SRA_id}_2_val_2.fq.gz")
    input:
        FWD="resources/atacseq/files/fastq/{SRA_id}/{SRA_id}_1.fastq.gz",
        REV="resources/atacseq/files/fastq/{SRA_id}/{SRA_id}_2.fastq.gz"
    log:
        "logs/trim_reads_atacseq/{SRA_id}.log"
    conda:
        "../../envs/trim_galore.yml"
    shell:
        """
        trim_galore --fastqc --gzip --paired \
                    -o resources/atacseq/files/fastq/{wildcards.SRA_id}/ \
                    {input.FWD} {input.REV} 2> {log}
        """

