rule get_COSMIC_mutations_coordinates:
    output:
        path_coding="resources/COSMIC/coordinates/coding/{gender_COSMIC}_{COSMIC_parameters}.bed",
        path_noncoding="resources/COSMIC/coordinates/noncoding/{gender_COSMIC}_{COSMIC_parameters}.bed"
    input:
        path_coding="resources/COSMIC/mutations/coding/{gender_COSMIC}_{COSMIC_parameters}.filtered.bed",
        path_noncoding="resources/COSMIC/mutations/noncoding/{gender_COSMIC}_{COSMIC_parameters}.filtered.bed"
    log:
        "logs/get_COSMIC_mutations_coordinates/{gender_COSMIC}_{COSMIC_parameters}.log"
    shell:
        """
        tail -n+2 {input.path_coding} | cut -f 7,15,16,17 | awk -F'\t' -v OFS='\t' '{{$2 = "chr"$2; print $2, $3, $4, $1}}' > {output.path_coding} 2>>{log}
        tail -n+2 {input.path_noncoding} | cut -f 8,11,12,13 | awk -F'\t' -v OFS='\t' '{{$2 = "chr"$2; print $2, $3, $4, $1}}' > {output.path_noncoding} 2>>{log}
        """