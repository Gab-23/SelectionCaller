rule sort_collapse_coding_mutations:
    output:
        "resources/COSMIC/mutations/coding/coding_mutations.filtered.unique.sorted.bed"
    input:
        "resources/COSMIC/mutations/coding/coding_mutations.filtered.bed"
    conda:
        "../envs/bedtools.yml"
    log:
        "logs/sort_collapse_coding_mutations/sort_collapse_coding_mutations.log"
    shell:
        """
	    bedtools sort -i {input} | uniq > {output} 2> {log}
        """