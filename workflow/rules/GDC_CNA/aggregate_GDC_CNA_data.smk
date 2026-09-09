rule aggregate_GDC_CNA_data:
    output:
        "resources/GDC/CNAs/aggregated/{tumor_type_GDC}_{gender_GDC}.bed"
    input:
        "resources/GDC/CNAs/checks/{tumor_type_GDC}_{gender_GDC}"
    log:
        "logs/process_GDC_CNA/{tumor_type_GDC}_{gender_GDC}.log"
    params:
        GDC_CNA_files_locations=config["GDC_CNA_files_locations"][0]
    shell:
        """
        FILES=$(ls {params.GDC_CNA_files_locations})
        for i in $FILES; do
        tail -n+2 $i | \
        cut -f 2,3,4,5,6 | \
        awk -F'\t' -v OFS='\t' '{{$1 = "chr"$1; print $1, $2, $3, $4, $5}}' >> {output} 2>> {log}; 
        done
        """