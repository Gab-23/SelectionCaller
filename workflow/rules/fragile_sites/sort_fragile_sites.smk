rule sort_fragile_sites:
    # NOTE (fix): fragile_sites_li_{assembly_code}.bed is a manually-curated source file (no
    # producing rule - see README) that was never run through bedtools sort at all. Every
    # other file compared against the backbone in this pipeline (backbone itself included,
    # via sort_backbone) is sorted with plain `bedtools sort` (no -g genome file), which uses
    # lexicographic chromosome order (chr1, chr10, chr11, ..., chr2, chr20, ...) rather than
    # natural numeric order - not itself a bug, just this pipeline's established convention,
    # but every file used together in the same bedtools call must share it. Confirmed directly
    # on a real cluster run: bedtools refused to proceed with "chromosome sort ordering...
    # inconsistent with other files" the moment it hit chr10 in the unsorted fragile sites file.
    output:
        "resources/fragile_sites/origin/fragile_sites_li_{assembly_code}.sorted.bed"
    input:
        "resources/fragile_sites/origin/fragile_sites_li_{assembly_code}.bed"
    log:
        "logs/sort_fragile_sites/{assembly_code}.log"
    conda:
        "../../envs/bedtools.yml"
    shell:
        """
        bedtools sort -i {input} > {output} 2> {log}
        """
