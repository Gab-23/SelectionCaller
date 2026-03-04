rule fastq_dump_atacseq:
    output:
        temp("resources/atacseq/files/fastq/{SRA_id}/{SRA_id}_1.fastq.gz"),
        temp("resources/atacseq/files/fastq/{SRA_id}/{SRA_id}_2.fastq.gz")
    input:
        "resources/atacseq/files/sra/{SRA_id}/{SRA_id}.sra"
    log:
        "logs/fastq_dump_atacseq/{SRA_id}.log"
    conda:
        "../envs/sra_toolkit.yml"
    shell:
        """
        fastq-dump --gzip \
                   --origfmt \
                   --split-files \
                   --skip-technical \
                   -O resources/atacseq/files/fastq/{wildcards.SRA_id}/ \
                   {input}
        """
