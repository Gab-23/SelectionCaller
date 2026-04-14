rule get_CNVkit_scaling_factors_atacseq:
    output:
        "resources/atacseq/files/cnvkit/{SRA_tumor_id}/{SRA_tumor_id}_{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}.scaling"
    input:
        path_backbone="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.filtered.sorted.bed",
        path_cnr="resources/atacseq/files/cnvkit/{SRA_tumor_id}/{SRA_tumor_id}_{assembly_code}_{mapq_thr}.processed"
    log:
        "logs/get_CNVkit_scaling_factors/{SRA_tumor_id}_{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}.log"
    conda:
        "../envs/bedtools_env.yml"
    shell:
        """
        bedtools closest -a {input.path_backbone} \
        -b <(bedtools sort -i {input.path_cnr}) \
        -wa -wb | \
        awk -v OFS='\t' '{{print $1,$2,$3,$4,$5"_"$6"_"$7,$8}}' \
        > {output}
        """
