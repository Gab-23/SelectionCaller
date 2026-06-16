rule prefetch_hic:
    output:
        temp("resources/hic/files/sra/{SRA_id}/{SRA_id}.sra")        
    log:
        "logs/prefetch_hic/{SRA_id}.log"
    shell:
        """
        mkdir -p resources/hic/files/sra/{wildcards.SRA_id}

        wget -O {output} \
            https://sra-pub-run-odp.s3.amazonaws.com/sra/{wildcards.SRA_id}/{wildcards.SRA_id} \
            2> {log}
        """