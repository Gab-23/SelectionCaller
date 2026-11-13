rule split_GDC_CNA_data:
    # NOTE (29/10/2026): divisione sul SEGNO della colonna Segment_Mean (quinta), non piu' con grep "-" sull'intera riga:
    # con l'identificativo del campione (UUID, contiene "-") nella sesta colonna, grep "-" metterebbe tutti i segmenti nelle
    # delezioni. Risultato identico a prima per i file a cinque colonne (dopo il filtro |Segment_Mean| > soglia nessun valore e' 0).
    output:
        path_ampl="resources/GDC/CNAs/ampl/{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}.bed",
        path_del="resources/GDC/CNAs/del/{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}.bed"
    input:
        "resources/GDC/CNAs/aggregated/{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}.filtered.sorted.bed"
    log:
        "logs/split_GDC_CNA_data/{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}.log"
    shell:
        """
        awk -F'\\t' '$5 > 0' {input} > {output.path_ampl} 2>{log}
        awk -F'\\t' '$5 < 0' {input} > {output.path_del} 2>>{log}
        """