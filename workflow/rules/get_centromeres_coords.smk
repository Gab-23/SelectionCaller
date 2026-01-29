rule get_centromeres_coords:
    output:
        "resources/reference_data/centromeres/{assembly_code}.bed"
    input:
        "resources/reference_data/cytobands/{assembly_code}_cytoBand.txt.gz"
    log:
        "logs/get_centromeres_coords/{assembly_code}.log"
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        grep acen {input} | bedtools merge > {output} 2> {log}
        """