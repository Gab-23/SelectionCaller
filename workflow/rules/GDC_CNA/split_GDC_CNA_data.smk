rule split_GDC_CNA_data:
    output:
        path_ampl="resources/GDC/CNAs/ampl/{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}.bed",
        path_del="resources/GDC/CNAs/del/{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}.bed"
    input:
        "resources/GDC/CNAs/aggregated/{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}.filtered.sorted.bed"
    log:
        "logs/split_GDC_CNA_data/{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}.log"
    shell:
        """
        grep -v "-" {input} > {output.path_ampl} 2>{log}
        grep "-" {input} > {output.path_del} 2>{log}
        """
