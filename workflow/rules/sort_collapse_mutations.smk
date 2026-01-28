rule sort_collapse_mutations:
    output:
        path_coding="resources/COSMIC/coordinates/coding/coordinates.sorted.bed",
        path_noncoding="resources/COSMIC/coordinates/noncoding/coordinates.sorted.bed"
    input:
        path_coding="resources/COSMIC/coordinates/coding/coordinates.bed",
        path_noncoding="resources/COSMIC/coordinates/noncoding/coordinates.bed"
    conda:
        "../envs/bedtools.yml"
    log:
        "logs/sort_collapse_mutations/sort_collapse_mutations.log"
    shell:
        """
	    sort -k 4 {input.path_coding} | uniq | bedtools sort > {output.path_coding} 2>> {log}
        sort -k 4 {input.path_noncoding} | uniq | bedtools sort > {output.path_noncoding} 2>> {log}
        """