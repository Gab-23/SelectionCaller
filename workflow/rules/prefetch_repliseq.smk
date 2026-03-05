rule prefetch_repliseq:
    output:
        temp("resources/repliseq/files/sra/{SRA_id}/{SRA_id}.sra")        
    log:
        "logs/prefetch_repliseq/{SRA_id}.log"
    shell:
        """
        mkdir -p resources/repliseq/files/sra/{wildcards.SRA_id}

        wget -O {output} \
            https://sra-pub-run-odp.s3.amazonaws.com/sra/{wildcards.SRA_id}/{wildcards.SRA_id} \
            2> {log}
        """