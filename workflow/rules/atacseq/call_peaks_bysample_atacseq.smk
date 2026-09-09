rule call_peaks_bysample_atacseq:
    output:
        pq_vals_path="resources/atacseq/files/peak_calling/QC/{SRA_id}/{assembly_code}_{mapq_thr}_{min_length_atacseq}_{q_val_thr_atacseq}_pq_vals",
        pileups_path="resources/atacseq/files/peak_calling/QC/{SRA_id}/{assembly_code}_{mapq_thr}_{min_length_atacseq}_{q_val_thr_atacseq}_pileups",
        analized_reads_path="resources/atacseq/files/peak_calling/QC/{SRA_id}/{assembly_code}_{mapq_thr}_{min_length_atacseq}_{q_val_thr_atacseq}_analyzed_reads",
        peaks_path="resources/atacseq/files/peak_calling/peaks/{SRA_id}/{assembly_code}_{mapq_thr}_{min_length_atacseq}_{q_val_thr_atacseq}.narrowPeak"
    input:
        "resources/atacseq/files/bam/{SRA_id}/{assembly_code}_{mapq_thr}.name.sorted.bam"
    log:
        "logs/call_peaks_conditionwise_atacseq/{SRA_id}_{assembly_code}_{mapq_thr}_{min_length_atacseq}_{q_val_thr_atacseq}.log"
    conda:
        "../../envs/call_peaks.yml"
    shell:
        """
        Genrich -j -D \
        -t {input} \
		-l {wildcards.min_length_atacseq} \
		-q {wildcards.q_val_thr_atacseq} \
		-f {output.pq_vals_path} \
		-k {output.pileups_path} \
		-b {output.analized_reads_path} \
		-o {output.peaks_path} 2> {log}
        """
