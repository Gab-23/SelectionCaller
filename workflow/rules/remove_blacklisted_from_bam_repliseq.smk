rule remove_blacklisted_from_bam_repliseq:
    output:
        "resources/repliseq/files/bam/{assembly_code}_{mapq_thr}_{SRA_id}.noblack.bam"
    input:
        bam="resources/repliseq/files/bam/{assembly_code}_{mapq_thr}_{SRA_id}.dedup.filtered.bam",
        blacklist="resources/blacklist/{assembly_code}.bed"
    log:
        "logs/remove_blacklisted_from_bam_repliseq/{assembly_code}_{mapq_thr}_{SRA_id}.log"
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        bedtools intersect -v -abam {input.bam} -b {input.blacklist} > {output}
        """