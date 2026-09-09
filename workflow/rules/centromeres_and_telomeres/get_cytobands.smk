rule get_cytobands:
    output:
        "resources/reference_data/cytobands/{assembly_code}_cytoBand.txt"
    log:
        "logs/get_cytobands/{assembly_code}.log"
    shell:
        """
        wget -qO- https://hgdownload.soe.ucsc.edu/goldenPath/{assembly_code}/database/cytoBand.txt.gz | gunzip | grep -v ["_MXY"] > {output} 2> {log}
        """