rule fix_mate_information_atacseq:
    output:
        "resources/atacseq/files/bam/{SRA_id}/{assembly_code}_{mapq_thr}.fixed.bam"
    input:
        "resources/atacseq/files/bam/{SRA_id}/{assembly_code}_{mapq_thr}.noblack.sorted.bam"
    log:
        "logs/fix_mate_information_atacseq/{assembly_code}_{mapq_thr}_{SRA_id}.log"
    conda:
        "../envs/picard.yml"
    shell:
        """
        picard FixMateInformation \
        I={input} \
        O={output} \
        ADD_MATE_CIGAR=true
        """