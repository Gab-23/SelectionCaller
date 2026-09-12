rule prefetch_hic:
    output:
        temp("resources/hic/files/sra/{SRA_id}/{SRA_id}.sra")        
    log:
        "logs/prefetch_hic/{SRA_id}.log"
    shell:
        # NOTE (fix): homologated with get_reference_genome's download-robustness fix -
        # --tries/--timeout for automatic retry on transient network issues, --continue so a
        # retry resumes instead of restarting from scratch.
        """
        mkdir -p resources/hic/files/sra/{wildcards.SRA_id}

        wget --tries=5 --timeout=60 --continue -O {output} \
            https://sra-pub-run-odp.s3.amazonaws.com/sra/{wildcards.SRA_id}/{wildcards.SRA_id} \
            2> {log}
        """
