rule prefetch_atacseq:
    output:
        "resources/atacseq/files/sra/{SRA_id}/{SRA_id}.sra"     
    log:
        "logs/prefetch_atacseq/{SRA_id}.log"
    conda:
        "../envs/sra_toolkit.yml"
    shell:
        """
        cd resources/atacseq/files/sra/
        prefetch {wildcards.SRA_id}
        """