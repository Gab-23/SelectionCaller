rule get_gencode_annotation:
    output:
        "resources/reference_data/gencode/gencode.v{gencode_version}.annotation.gff3.gz",
        "resources/reference_data/gencode/gencode.v{gencode_version}.annotation.gtf.gz"
    log:
        "logs/get_gencode_annotation/v{gencode_version}.log"
    shell:
        """
        wget -qP resources/reference_data/gencode/ https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{wildcards.gencode_version}/gencode.v{wildcards.gencode_version}.annotation.gff3.gz 2>> {log}
        wget -qP resources/reference_data/gencode/ https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{wildcards.gencode_version}/gencode.v{wildcards.gencode_version}.annotation.gtf.gz 2>> {log}
        """