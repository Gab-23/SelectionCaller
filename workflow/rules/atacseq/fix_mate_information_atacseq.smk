rule fix_mate_information_atacseq:
    output:
        # NOTE (fix, corrected direction): the filename genuinely needs the {SRA_id} prefix -
        # cnvkit.py batch derives its internal sample ID from each bam's basename alone
        # (ignoring the directory), so multiple samples sharing an identical basename (as
        # after an earlier, well-intentioned but wrong fix that dropped this prefix) makes
        # CNVkit reject the whole batch with "Duplicate sample ID", confirmed with a real
        # Snakemake run. The 3 other consumers of this file (name_sort_final_bam_atacseq,
        # compute_insert_size_metrics_atacseq, compute_coverage_atacseq) only ever process one
        # sample's bam at a time, so they don't care either way - updated to match instead.
        "resources/atacseq/files/bam/{SRA_id}/{SRA_id}_{assembly_code}_{mapq_thr}.fixed.bam"
    input:
        "resources/atacseq/files/bam/{SRA_id}/{assembly_code}_{mapq_thr}.noblack.sorted.bam"
    log:
        "logs/fix_mate_information_atacseq/{assembly_code}_{mapq_thr}_{SRA_id}.log"
    conda:
        "../../envs/picard.yml"
    shell:
        """
        picard FixMateInformation \
        I={input} \
        O={output} \
        ADD_MATE_CIGAR=true 2>{log}
        """
