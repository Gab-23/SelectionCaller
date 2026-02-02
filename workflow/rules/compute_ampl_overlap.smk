rule compute_ampl_overlap:
    output:
        "resources/num_overlaps/CNAs/ampl/{assembly_code}_{window_mode}_{window_size}_{step}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}"
    input:
        path_CNA="resources/GDC/CNAs/ampl/{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}.bed",
        path_backbone="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.filtered.sorted.bed"
    log:
        "logs/compute_GDcompute_ampl_overlapC_CNA_overlap/{assembly_code}_{window_mode}_{window_size}_{step}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}.log"
    conda:
        "../envs/bedtools.yml"
    shell:
        """
        bedtools intersect -a {input.path_backbone} -b {input.path_CNA} -c > {output} 2>> {log}
        """