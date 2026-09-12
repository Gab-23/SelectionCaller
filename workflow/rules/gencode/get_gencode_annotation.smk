rule get_gencode_annotation:
    output:
        "resources/reference_data/gencode/gencode.v{gencode_version}.annotation.gff3.gz",
        "resources/reference_data/gencode/gencode.v{gencode_version}.annotation.gtf.gz"
    log:
        "logs/get_gencode_annotation/v{gencode_version}.log"
    shell:
        # NOTE (fix): homologated with get_reference_genome's download-robustness fix -
        # dropped -q (was hiding wget's own diagnostics), added --tries/--timeout for
        # automatic retry, --continue so a retry resumes instead of restarting.
        """
        wget --tries=5 --timeout=60 --continue -P resources/reference_data/gencode/ https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{wildcards.gencode_version}/gencode.v{wildcards.gencode_version}.annotation.gff3.gz 2>> {log}
        wget --tries=5 --timeout=60 --continue -P resources/reference_data/gencode/ https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{wildcards.gencode_version}/gencode.v{wildcards.gencode_version}.annotation.gtf.gz 2>> {log}
        """
