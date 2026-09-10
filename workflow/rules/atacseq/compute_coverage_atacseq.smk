rule compute_coverage_atacseq:
    output:
        "resources/atacseq/files/coverage/{SRA_id}/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}.bed"
    input:
        backbone="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.mask.filtered.sorted.bed",
        bam="resources/atacseq/files/bam/{SRA_id}/{SRA_id}_{assembly_code}_{mapq_thr}.fixed.bam"
    log:
        "logs/compute_coverage_atacseq/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{SRA_id}.log"
    conda:
        "../../envs/bedtools.yml"
    shell:
        """
        bedtools coverage -counts -a {input.backbone} -b {input.bam} > {output}
        """
