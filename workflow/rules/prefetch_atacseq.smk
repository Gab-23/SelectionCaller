rule prefetch_atacseq:
    output:
        temp("resources/atacseq/files/sra/{SRA_id}/{SRA_id}.sra")
    log:
        "logs/prefetch_atacseq/{SRA_id}.log"
    shell:
        """
        mkdir -p resources/atacseq/files/sra/{wildcards.SRA_id}

        wget -O {output} \
            https://sra-pub-run-odp.s3.amazonaws.com/sra/{wildcards.SRA_id}/{wildcards.SRA_id} \
            2> {log}
        """