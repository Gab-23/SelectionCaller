rule compute_del_overlap:
    # NOTE (29/10/2026): numero di PAZIENTI DISTINTI con un segmento che tocca il bin (bedtools map, count_distinct sulla sesta
    # colonna = campione), non piu' numero di SEGMENTI (bedtools intersect -c): un paziente con un punto di rottura dentro il bin
    # contava piu' volte. Stesso formato di uscita (backbone + conteggio); process_GDC_CNA_overlap divide per il numero di pazienti.
    # Richiede che backbone e segmenti siano ordinati allo stesso modo (entrambi con bedtools sort, come gia' nella pipeline).
    output:
        "resources/num_overlaps/CNAs/del/{assembly_code}_{window_mode}_{window_size}_{step}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}"
    input:
        path_CNA="resources/GDC/CNAs/del/{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}.bed",
        path_backbone="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.mask.filtered.sorted.bed"
    log:
        "logs/compute_del_overlap/{assembly_code}_{window_mode}_{window_size}_{step}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}.log"
    conda:
        "../../envs/bedtools.yml"
    shell:
        """
        bedtools map -a {input.path_backbone} -b {input.path_CNA} -c 6 -o count_distinct -null 0 > {output} 2>> {log}
        """