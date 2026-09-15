rule sort_raw_backbone_hic:
    # NOTE (revised Option B design): instead of generating a separate continuous bin table
    # via `cooler makebins` (a parallel mechanism to the rest of the pipeline's own backbone
    # construction), this reuses the RAW, pre-filter backbone
    # ("resources/backbone/{...}.bed", from get_backbone_file.R) - which is already a
    # continuous, gap-free tiling of the whole genome (built directly from
    # {assembly_code}.clean.chrom.sizes, before filter_backbone ever removes the blacklist)
    # and already carries the same binID convention as the rest of the pipeline. This rule
    # only sorts it (same bedtools sort, no -g, as sort_backbone.smk already does for the
    # filtered version) - a small, dedicated rule rather than modifying sort_backbone.smk
    # itself, to avoid any risk to its many other existing consumers.
    output:
        "resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.sorted.bed"
    input:
        "resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.bed"
    log:
        "logs/sort_raw_backbone_hic/{assembly_code}_{window_mode}_{window_size}_{step}.log"
    conda:
        "../../envs/bedtools.yml"
    shell:
        """
        bedtools sort -i {input} > {output} 2> {log}
        """