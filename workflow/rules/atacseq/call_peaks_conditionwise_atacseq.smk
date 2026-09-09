rule call_peaks_conditionwise_atacseq:
    output:
        pq_vals_path="resources/atacseq/files/peaks/QC/{condition}/{assembly_code}_{mapq_thr}_{min_length_atacseq}_{q_val_thr_atacseq}_pq_vals,
        pileups_path="resources/atacseq/files/peaks/QC/{condition}/{assembly_code}_{mapq_thr}_{min_length_atacseq}_{q_val_thr_atacseq}_pileups,
        analized_reads_path="resources/atacseq/files/peaks/QC/{condition}/{assembly_code}_{mapq_thr}_{min_length_atacseq}_{q_val_thr_atacseq}_analyzed_reads,
        peaks_path="resources/atacseq/files/peaks/{condition}/{assembly_code}_{mapq_thr}_{min_length_atacseq}_{q_val_thr_atacseq}.narrowPeak
    input:
        check="resources/atacseq/checks/{assembly_code}_{mapq_thr}_name_sorting",
        metadata=atacseq_selected_samples_location
    log:
        "logs/call_peaks_conditionwise_atacseq/{condition}_{assembly_code}_{mapq_thr}_{min_length_atacseq}_{q_val_thr_atacseq}.log"
    conda:
        "../../envs/call_peaks.yml"
    shell:
        """
        GENRICH_INPUT=$(awk -v condition="{wildcards.condition}" '$6 == condition {print $1}' {input.metadata} | \
                        while read id; do
                            find resources/atacseq/files/bam -type f -path "*/${id}/*.name.sorted.bam"
                        done | \
                        tr '\n' ' ')

        Genrich -j -D \
        -t $GENRICH_INPUT \
		-l {wildcards.min_length_atacseq} \
		-q {wildcards.q_val_thr_atacseq} \
		-f {output.pq_vals_path} \
		-k {output.pileups_path} \
		-b {output.analized_reads_path} \
		-o {output.peaks_path} 2> {log}
        """