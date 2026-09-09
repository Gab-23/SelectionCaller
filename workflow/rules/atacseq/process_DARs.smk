rule process_DARs:
    output:
        "data/atacseq/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}"
    input:
        "resources/atacseq/files/DARs/woCNAs/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}"
    log:
        "logs/process_DARs/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}.log"
    shell:
        """
        awk -F '\t' -v OFS='\t' '{print $6,$1}' {input} > {output} 2> {log}
        """