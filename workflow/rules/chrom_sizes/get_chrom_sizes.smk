rule get_chrom_sizes:
    output:
        "resources/reference_data/chrom_sizes/{assembly_code}.chrom.sizes"
    log:
        "logs/get_chrom_sizes/{assembly_code}.log"
    shell:
        # NOTE (fix): homologated with get_reference_genome's download-robustness fix -
        # dropped -q (was hiding wget's own diagnostics), added --tries/--timeout for
        # automatic retry, --continue so a retry resumes instead of restarting.
        """
        wget --tries=5 --timeout=60 --continue -P resources/reference_data/chrom_sizes/ https://hgdownload.soe.ucsc.edu/goldenpath/{wildcards.assembly_code}/bigZips/{wildcards.assembly_code}.chrom.sizes 2> {log}
        """
