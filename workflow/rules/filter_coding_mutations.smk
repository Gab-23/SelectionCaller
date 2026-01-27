rule filter_coding_mutations:
    output:
        "resources/COSMIC/mutations/coding/coding_mutations.filtered.bed"
    input:
        path_sample_ids="resources/COSMIC/samples/final_samples.tsv",
        path_mutation="resources/COSMIC/mutations/coding/Cosmic_GenomeScreensMutant_Tsv_v103_GRCh38/Cosmic_GenomeScreensMutant_v103_GRCh38.tsv"
    log:
        "logs/map_COSMIC_sample_mutation_IDs/map_COSMIC_sample_mutation_IDs.log"
    shell:
        """
        grep -Fw -f {input.path_sample_ids} {input.path_mutation} | cut -f 15,16,17,12 | awk -F'\t' -v OFS='\t' '{{$2 = "chr"$2; print $2, $3, $4, $1}}' >> {output} 2>>{log}
        """