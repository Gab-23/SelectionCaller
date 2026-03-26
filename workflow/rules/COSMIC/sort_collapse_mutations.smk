rule sort_collapse_mutations:
    output:
        path_coding="resources/COSMIC/coordinates/coding/{gender_COSMIC}_{COSMIC_parameters}.sorted.bed",
        path_noncoding="resources/COSMIC/coordinates/noncoding/{gender_COSMIC}_{COSMIC_parameters}.sorted.bed"
    input:
        path_coding="resources/COSMIC/coordinates/coding/{gender_COSMIC}_{COSMIC_parameters}.bed",
        path_noncoding="resources/COSMIC/coordinates/noncoding/{gender_COSMIC}_{COSMIC_parameters}.bed"
    conda:
        "../envs/bedtools.yml"
    log:
        "logs/sort_collapse_mutations/{gender_COSMIC}_{COSMIC_parameters}.log"
    shell:
        """
	    sort -k 4 {input.path_coding} | uniq | bedtools sort > {output.path_coding} 2>> {log}
        sort -k 4 {input.path_noncoding} | uniq | bedtools sort > {output.path_noncoding} 2>> {log}
        """