rule sort_GDC_CNA_data:
    output:
        "resources/GDC/CNAs/aggregated/{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}.filtered.sorted.bed"
    input:
        "resources/GDC/CNAs/aggregated/{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}.filtered.bed"
    log:
        "logs/sort_GDC_CNA_data/{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}.log"
    conda:
        "../../envs/bedtools.yml"
    shell:
        """
        bedtools sort -i {input} > {output} 2> {log}
        """