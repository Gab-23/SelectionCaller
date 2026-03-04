rule prefetch_repliseq:
    output:
        "resources/repliseq/files/sra/{SRA_id}/{SRA_id}.sra"        
    log:
        "logs/prefetch_repliseq/{SRA_id}.log"
    conda:
        "../envs/sra_toolkit.yml"
    shell:
        """
        cd resources/repliseq/files/sra/
        prefetch --force all --transport https {wildcards.SRA_id}
        """