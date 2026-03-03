rule compute_insert_size_metrics_atacseq:
    output:
        metrics="results/atacseq/IS/{SRA_id}/{assembly_code}_{mapq_thr}",
        hist="results/atacseq/IS/{SRA_id}/{assembly_code}_{mapq_thr}.pdf"
    input:
        "resources/atacseq/files/bam/{SRA_id}/{assembly_code}_{mapq_thr}.fixed.bam"
    log:
        "logs/compute_insert_size_metrics_atacseq/{assembly_code}_{mapq_thr}_{SRA_id}.log"
    conda:
        "../envs/picard.yml"
    shell:
        """
        picard CollectInsertSizeMetrics \
        I={input} \
        O={output.metrics} \
        H={output.hist}
        """