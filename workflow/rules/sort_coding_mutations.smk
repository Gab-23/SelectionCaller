rule sort_coding_mutations:
    output:
        "resources/COSMIC/mutations/coding/mutations_w_description.sorted.bed"
    input:
        "resources/COSMIC/mutations/coding/mutations_w_description.bed"
    conda:
        "../envs/bedtools.yml"
    log:
        "logs/sort_coding_mutations/sort_coding_mutations.log"
    shell:
        """
	    bedtools sort -i {input} > {output} 2> {log}
        """