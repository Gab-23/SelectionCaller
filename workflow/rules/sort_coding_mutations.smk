rule sort_coding_mutations:
    output:
        "resources/COSMIC/mutations/coding/{gender_COSMIC}_{COSMIC_parameters}_w_description.sorted.bed"
    input:
        "resources/COSMIC/mutations/coding/{gender_COSMIC}_{COSMIC_parameters}_w_description.bed"
    conda:
        "../envs/bedtools.yml"
    log:
        "logs/sort_coding_mutations/{gender_COSMIC}_{COSMIC_parameters}.log"
    shell:
        """
	    bedtools sort -i {input} > {output} 2> {log}
        """