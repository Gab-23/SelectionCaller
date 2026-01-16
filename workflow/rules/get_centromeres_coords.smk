rule get_centromeres_coords:
    output:
        "resources/reference_data/centromeres/{assembly_code}.bed"
    input:
        "resources/reference_data/cytobands/{assembly_code}_cytoBand.txt.gz"
    conda:
        "../envs/bedtools.yml"
    log:
        "logs/get_centromeres_coords/{assembly_code}.log"
    shell:
        """
        grep acen {input} | bedtools merge > {output} 2> {log}
        """