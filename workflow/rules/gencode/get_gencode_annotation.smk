rule get_gencode_annotation:
    # GENCODE gene annotation (GFF3 + GTF).
    #   hg38: primary release, release_{N}/gencode.v{N}.annotation.*
    #   hg19: the same GRCh38 release back-mapped to GRCh37 with gencode-backmap,
    #         release_{N}/GRCh37_mapping/gencode.v{N}lift37.annotation.*
    # The snakefile turns gencode_version into "{N}lift37" when assembly_code == "hg19", so
    # the wildcard IS GENCODE's own file version and every downstream path / dataset name
    # records which annotation was used (hg19 and hg38 files never collide).
    # lift37 notes (GENCODE GRCh37_mapping README): IDs of back-mapped features carry a
    # mapping version suffix (e.g. ENSG00000223972.5_1); genes that could not be mapped
    # consistently are taken from v19 (possibly with a different gene_id) or are missing;
    # some genes are placed on unplaced scaffolds (dropped by the backbone intersection).
    output:
        "resources/reference_data/gencode/gencode.v{gencode_version}.annotation.gff3.gz",
        "resources/reference_data/gencode/gencode.v{gencode_version}.annotation.gtf.gz"
    params:
        url_dir=lambda wildcards: (
            "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_"
            + wildcards.gencode_version.replace("lift37", "")
            + ("/GRCh37_mapping" if wildcards.gencode_version.endswith("lift37") else "")
        )
    log:
        "logs/get_gencode_annotation/v{gencode_version}.log"
    shell:
        """
        wget --tries=5 --timeout=60 --continue -P resources/reference_data/gencode/ {params.url_dir}/gencode.v{wildcards.gencode_version}.annotation.gff3.gz 2>> {log}
        wget --tries=5 --timeout=60 --continue -P resources/reference_data/gencode/ {params.url_dir}/gencode.v{wildcards.gencode_version}.annotation.gtf.gz 2>> {log}
        """