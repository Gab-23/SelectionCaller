rule compute_coverage_repliseq:
    output:
        "resources/repliseq/files/coverage/{SRA_id}/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}.bed"
    input:
        backbone="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.filtered.sorted.bed",
        bam="resources/repliseq/files/bam/{SRA_id}/{assembly_code}_{mapq_thr}.noblack.sorted.bam"
    log:
        "logs/compute_coverage_repliseq/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{SRA_id}.log"
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        bedtools coverage -counts -a {input.backbone} -b {input.bam} > {output}
        """
