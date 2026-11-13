rule aggregate_GDC_CNA_data:
    # NOTE (29/10/2026): si aggregano SOLO i file elencati nel manifest della coorte (resources/GDC/CNAs/files/<id>/<filename>),
    # non piu' tutti i file della cartella (config GDC_CNA_files_locations): con piu' coorti scaricate nella stessa cartella,
    # l'aggregato di una coorte conteneva anche i segmenti delle altre (LUAD_ALL: 1572 file invece di 497, frequenze fino a 2,5).
    # Se un file del manifest manca, la regola fallisce.
    # NOTE (29/10/2026, frequenze per PAZIENTE): si tiene anche l'identificativo del campione (colonna GDC_Aliquot, la prima del
    # file .seg), come SESTA colonna: chr, start, end, Num_Probes, Segment_Mean, campione. Serve a contare i pazienti distinti per
    # bin in compute_ampl_overlap / compute_del_overlap. Le prime cinque colonne sono invariate.
    output:
        "resources/GDC/CNAs/aggregated/{tumor_type_GDC}_{gender_GDC}.bed"
    input:
        check="resources/GDC/CNAs/checks/{tumor_type_GDC}_{gender_GDC}",
        manifest="resources/GDC/CNAs/{tumor_type_GDC}_{gender_GDC}_manifest"
    log:
        "logs/process_GDC_CNA/{tumor_type_GDC}_{gender_GDC}.log"
    shell:
        """
        : > {output}
        n=0
        while IFS=$'\\t' read -r id fname rest; do
            f="resources/GDC/CNAs/files/$id/$fname"
            if [ ! -f "$f" ]; then echo "file del manifest mancante: $f" >> {log}; exit 1; fi
            tail -n+2 "$f" | \
            cut -f 1,2,3,4,5,6 | \
            awk -F'\\t' -v OFS='\\t' '{{print "chr"$2, $3, $4, $5, $6, $1}}' >> {output} 2>> {log}
            n=$((n + 1))
        done < <(tail -n+2 {input.manifest})
        echo "file aggregati dal manifest: $n" >> {log}
        """