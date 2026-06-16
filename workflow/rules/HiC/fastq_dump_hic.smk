rule fastq_dump_hic:
    output:
        temp("resources/hic/files/fastq/{SRA_id}/{SRA_id}_1.fastq.gz"),
        temp("resources/hic/files/fastq/{SRA_id}/{SRA_id}_2.fastq.gz")
    input:
        "resources/hic/files/sra/{SRA_id}/{SRA_id}.sra"
    log:
        "logs/fastq_dump_hic/{SRA_id}.log"
    conda:
        "../envs/sra_toolkit.yml"
    shell:
        """
        fastq-dump --gzip \
                   --origfmt \
                   --split-files \
                   --skip-technical \
                   -O resources/hic/files/fastq/{wildcards.SRA_id}/ \
                   {input}
        """
