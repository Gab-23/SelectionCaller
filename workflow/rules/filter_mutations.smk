rule filter_mutations:
    output:
        path_coding="resources/COSMIC/mutations/coding/mutations.filtered.bed",
        path_noncoding="resources/COSMIC/mutations/noncoding/mutations.filtered.bed"
    input:
        path_sample_ids="resources/COSMIC/samples/final_samples.tsv",
        path_coding="resources/COSMIC/mutations/coding/Cosmic_GenomeScreensMutant_Tsv_v103_GRCh38/Cosmic_GenomeScreensMutant_v103_GRCh38.tsv",
        path_noncoding="resources/COSMIC/mutations/noncoding/Cosmic_NonCodingVariants_Tsv_v103_GRCh38/Cosmic_NonCodingVariants_v103_GRCh38.tsv"
    log:
        "logs/filter_mutations/filter_mutations.log"
    shell:
        """
        head -n 1 {input.path_coding} >> {output.path_coding} 2>> {log}
        head -n 1 {input.path_noncoding} >> {output.path_noncoding} 2>> {log}
        grep -Fw -f {input.path_sample_ids} {input.path_coding} >> {output.path_coding} 2>> {log}
        grep -Fw -f {input.path_sample_ids} {input.path_noncoding} >> {output.path_noncoding} 2>> {log}
        """