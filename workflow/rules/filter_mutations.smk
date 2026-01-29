rule filter_mutations:
    output:
        path_coding="resources/COSMIC/mutations/coding/{gender_COSMIC}_{COSMIC_parameters}.filtered.bed",
        path_noncoding="resources/COSMIC/mutations/noncoding/{gender_COSMIC}_{COSMIC_parameters}.filtered.bed"
    input:
        path_sample_ids="resources/COSMIC/samples/{gender_COSMIC}_{COSMIC_parameters}",
        path_coding="resources/COSMIC/mutations/coding/Cosmic_GenomeScreensMutant_Tsv_v103_GRCh38/Cosmic_GenomeScreensMutant_v103_GRCh38.tsv",
        path_noncoding="resources/COSMIC/mutations/noncoding/Cosmic_NonCodingVariants_Tsv_v103_GRCh38/Cosmic_NonCodingVariants_v103_GRCh38.tsv"
    log:
        "logs/filter_mutations/{gender_COSMIC}_{COSMIC_parameters}.log"
    shell:
        """
        head -n 1 {input.path_coding} >> {output.path_coding} 2>> {log}
        head -n 1 {input.path_noncoding} >> {output.path_noncoding} 2>> {log}
        grep -Fw -f {input.path_sample_ids} {input.path_coding} >> {output.path_coding} 2>> {log}
        grep -Fw -f {input.path_sample_ids} {input.path_noncoding} >> {output.path_noncoding} 2>> {log}
        """