rule fastq_dump_repliseq:
    output:
        temp("resources/repliseq/files/fastq/{SRA_id}/{SRA_id}_1.fastq.gz"),
        temp("resources/repliseq/files/fastq/{SRA_id}/{SRA_id}_2.fastq.gz")
    input:
        "resources/repliseq/files/sra/{SRA_id}/{SRA_id}.sra"
    log:
        "logs/fastq_dump_repliseq/{SRA_id}.log"
    conda:
        "../envs/sra_toolkit.yml"
    shell:
        """
        fastq-dump --gzip \
                   --origfmt \
                   --split-files \
                   --skip-technical \
                   -O resources/repliseq/files/fastq/{wildcards.SRA_id}/ \
                   {input}
        """
