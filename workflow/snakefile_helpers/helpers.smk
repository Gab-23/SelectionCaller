from snakemake.io import expand

rule all_prefetch_repliseq:
    input:
        get_prefetch_outputs_repliseq
    output:
        all_done="resources/repliseq/checks/prefetch"
    shell:
        """
        echo {input} >> {output}
        """

rule all_fastqdump_repliseq:
    input:
        get_fastqdump_outputs_repliseq
    output:
        all_done="resources/repliseq/checks/fastq_dump"
    shell:
        """
        echo {input} >> {output}
        """

rule all_trimgalore_repliseq:
    input:
        get_trim_galore_outputs_repliseq
    output:
        all_done="resources/repliseq/checks/trim_galore"
    shell:
        """
        echo {input} >> {output}
        """

rule all_alignment_repliseq:
    input:
        get_alignment_outputs_repliseq
    output:
        all_done="resources/repliseq/checks/{assembly_code}_alignment"
    shell:
        """
        echo {input} >> {output}
        """

rule all_alignment_sorting_repliseq:
    input:
        get_alignment_sorting_outputs_repliseq
    output:
        all_done="resources/repliseq/checks/{assembly_code}_alignment_sorting"
    shell:
        """
        echo {input} >> {output}
        """

rule all_alignment_indexing_repliseq:
    input:
        get_alignment_indexing_outputs_repliseq
    output:
        all_done="resources/repliseq/checks/{assembly_code}_alignment_indexing"
    shell:
        """
        echo {input} >> {output}
        """

rule all_filtering_repliseq:
    input:
        get_filtering_outputs_repliseq
    output:
        all_done="resources/repliseq/checks/{assembly_code}_alignment_filtering"
    shell:
        """
        echo {input} >> {output}
        """

rule all_sorting_filtered_alignment_repliseq:
    input:
        get_sorting_filtered_alignment_outputs_repliseq
    output:
        all_done="resources/repliseq/checks/{assembly_code}_filtered_alignment_sorting"
    shell:
        """
        echo {input} >> {output}
        """

rule all_indexing_filtered_alignment_repliseq:
    input:
        get_indexing_filtered_alignment_outputs_repliseq
    output:
        all_done="resources/repliseq/checks/{assembly_code}_filtered_alignment_indexing"
    shell:
        """
        echo {input} >> {output}
        """

rule all_bamqc_repliseq:
    input:
        get_bamqc_outputs_repliseq
    output:
        all_done="resources/repliseq/checks/{assembly_code}_bamqc"
    shell:
        """
        echo {input} >> {output}
        """

rule all_deduplicated_alignment_repliseq:
    input:
        get_deduplicated_alignment_outputs_repliseq
    output:
        all_done="resources/repliseq/checks/{assembly_code}_alignment_deduplication"
    shell:
        """
        echo {input} >> {output}
        """

rule all_filtered_deduplicated_alignment_repliseq:
    input:
        get_filtered_deduplicated_alignment_outputs_repliseq
    output:
        all_done="resources/repliseq/checks/{assembly_code}_{mapq_thr}_deduplicated_alignment_filtering"
    shell:
        """
        echo {input} >> {output}
        """

rule all_remove_blacklisted_alignment_repliseq:
    input:
        get_remove_blacklisted_alignment_outputs_repliseq
    output:
        all_done="resources/repliseq/checks/{assembly_code}_{mapq_thr}_noblacklist_alignment"
    shell:
        """
        echo {input} >> {output}
        """

rule all_sort_blacklisted_alignment_repliseq:
    input:
        get_sort_noblacklisted_alignment_outputs_repliseq
    output:
        all_done="resources/repliseq/checks/{assembly_code}_{mapq_thr}_noblacklist_alignment_sorting"
    shell:
        """
        echo {input} >> {output}
        """

rule all_index_blacklisted_alignment_repliseq:
    input:
        get_index_noblacklisted_alignment_outputs_repliseq
    output:
        all_done="resources/repliseq/checks/{assembly_code}_{mapq_thr}_noblacklist_alignment_indexing"
    shell:
        """
        echo {input} >> {output}
        """

rule all_get_coverage_repliseq:
    input:
        get_coverage_outputs_repliseq
    output:
        all_done="resources/repliseq/checks/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_coverage"
    shell:
        """
        echo {input} >> {output}
        """

rule all_filter_coverage_repliseq:
    input:
        filter_coverage_outputs_repliseq
    output:
        all_done="resources/repliseq/checks/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_coverage_filtering"
    shell:
        """
        echo {input} >> {output}
        """

rule all_prefetch_atacseq:
    input:
        get_prefetch_outputs_atacseq
    output:
        all_done="resources/atacseq/checks/prefetch"
    shell:
        """
        echo {input} >> {output}
        """

rule all_fastqdump_atacseq:
    input:
        get_fastqdump_outputs_atacseq
    output:
        all_done="resources/atacseq/checks/fastq_dump"
    shell:
        """
        echo {input} >> {output}
        """

rule all_trimgalore_atacseq:
    input:
        get_trim_galore_outputs_atacseq
    output:
        all_done="resources/atacseq/checks/trim_galore"
    shell:
        """
        echo {input} >> {output}
        """

rule all_alignment_atacseq:
    input:
        get_alignment_outputs_atacseq
    output:
        all_done="resources/atacseq/checks/{assembly_code}_alignment"
    shell:
        """
        echo {input} >> {output}
        """

rule all_alignment_sorting_atacseq:
    input:
        get_alignment_sorting_outputs_atacseq
    output:
        all_done="resources/atacseq/checks/{assembly_code}_alignment_sorting"
    shell:
        """
        echo {input} >> {output}
        """

rule all_alignment_indexing_atacseq:
    input:
        get_alignment_indexing_outputs_atacseq
    output:
        all_done="resources/atacseq/checks/{assembly_code}_alignment_indexing"
    shell:
        """
        echo {input} >> {output}
        """

rule all_filtering_atacseq:
    input:
        get_filtering_outputs_atacseq
    output:
        all_done="resources/atacseq/checks/{assembly_code}_alignment_filtering"
    shell:
        """
        echo {input} >> {output}
        """

rule all_sorting_filtered_alignment_atacseq:
    input:
        get_sorting_filtered_alignment_outputs_atacseq
    output:
        all_done="resources/atacseq/checks/{assembly_code}_filtered_alignment_sorting"
    shell:
        """
        echo {input} >> {output}
        """

rule all_indexing_filtered_alignment_atacseq:
    input:
        get_indexing_filtered_alignment_outputs_atacseq
    output:
        all_done="resources/atacseq/checks/{assembly_code}_filtered_alignment_indexing"
    shell:
        """
        echo {input} >> {output}
        """

rule all_bamqc_atacseq:
    input:
        get_bamqc_outputs_atacseq
    output:
        all_done="resources/atacseq/checks/{assembly_code}_bamqc"
    shell:
        """
        echo {input} >> {output}
        """

rule all_deduplicated_alignment_atacseq:
    input:
        get_deduplicated_alignment_outputs_atacseq
    output:
        all_done="resources/atacseq/checks/{assembly_code}_alignment_deduplication"
    shell:
        """
        echo {input} >> {output}
        """

rule all_filtered_deduplicated_alignment_atacseq:
    input:
        get_filtered_deduplicated_alignment_outputs_atacseq
    output:
        all_done="resources/atacseq/checks/{assembly_code}_{mapq_thr}_deduplicated_alignment_filtering"
    shell:
        """
        echo {input} >> {output}
        """

rule all_remove_blacklisted_alignment_atacseq:
    input:
        get_remove_blacklisted_alignment_outputs_atacseq
    output:
        all_done="resources/atacseq/checks/{assembly_code}_{mapq_thr}_noblacklist_alignment"
    shell:
        """
        echo {input} >> {output}
        """

rule all_sort_blacklisted_alignment_atacseq:
    input:
        get_sort_noblacklisted_alignment_outputs_atacseq
    output:
        all_done="resources/atacseq/checks/{assembly_code}_{mapq_thr}_noblacklist_alignment_sorting"
    shell:
        """
        echo {input} >> {output}
        """

rule all_index_blacklisted_alignment_atacseq:
    input:
        get_index_noblacklisted_alignment_outputs_atacseq
    output:
        all_done="resources/atacseq/checks/{assembly_code}_{mapq_thr}_noblacklist_alignment_indexing"
    shell:
        """
        echo {input} >> {output}
        """

rule all_fix_mate_information_atacseq:
    input:
        get_fixed_alignment_outputs_atacseq
    output:
        all_done="resources/atacseq/checks/{assembly_code}_{mapq_thr}_mate_information_fixing"
    shell:
        """
        echo {input} >> {output}
        """

rule all_get_is_metrics_atacseq:
    input:
        get_is_metrics_outputs_atacseq
    output:
        all_done="results/atacseq/checks/{assembly_code}_{mapq_thr}_is_metrics_getting"
    shell:
        """
        echo {input} >> {output}
        """

rule all_get_coverage_atacseq:
    input:
        get_coverage_outputs_atacseq
    output:
        all_done="resources/atacseq/checks/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_coverage"
    shell:
        """
        echo {input} >> {output}
        """

rule all_get_CNVkit_outputs_atacseq:
    input:
        get_CNVkit_outputs_atacseq
    output:
        all_done="resources/atacseq/checks/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_CNVkit"
    shell:
        """
	echo {input} >> {output}
	"""
rule all_process_CNVkit_outputs_atacseq:
    input:
        process_CNVkit_outputs_atacseq
    output:
        all_done="resources/atacseq/checks/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_processed_CNVkit"
    shell:
        """
        echo {input} >> {output}
        """

rule all_get_CNVkit_scaling_factors_atacseq:
    input:
        get_CNVkit_scaling_factors_atacseq
    output:
        all_done="resources/atacseq/checks/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_scaling_factors"
    shell:
        """
        echo {input} >> {output}
        """


#rule all_name_sort_final_bams_atacseq:
#    input:
#        get_final_bam_name_sorted_outputs_atacseq
#    output:
#        all_done="resources/atacseq/checks/{assembly_code}_{mapq_thr}_name_sorting"
#    shell:
#        """
#        echo {input} >> {output}
#        """
#
#rule all_peak_calling_atacseq:
#    input:
#        get_peak_calling_outputs_atacseq
#    output:
#        all_done="resources/atacseq/checks/{assembly_code}_{mapq_thr}_{min_length_atacseq}_{q_val_thr_atacseq}_peak_calling"
#    shell:
#        """
#        echo {input} >> {output}
#        """

#rule all_peak_processing_atacseq:
#    input:
#        get_peak_processing_outputs_atacseq
#    output:
#        all_done="resources/atacseq/checks/{assembly_code}_{mapq_thr}_{min_length_atacseq}_{q_val_thr_atacseq}_peak_processing"
#    shell:
#        """
#        echo {input} >> {output}
#        """
