rule explore_cohort_atacseq:
    output:
        path_samples_barplot="results/plots/atacseq/exploratory/samples_barplot/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{oncotree_code_DepMap}_{gender_DepMap}_{tumor_type_GDC}_{gender_GDC}_{gender_COSMIC}_{COSMIC_parameters}.pdf",
        path_plot_samples_per_study_by_condition="results/plots/atacseq/exploratory/samples_per_study_by_condition/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{oncotree_code_DepMap}_{gender_DepMap}_{tumor_type_GDC}_{gender_GDC}_{gender_COSMIC}_{COSMIC_parameters}.pdf",
        path_plot_subtype_distribution="results/plots/atacseq/exploratory/subtype_distribution/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{oncotree_code_DepMap}_{gender_DepMap}_{tumor_type_GDC}_{gender_GDC}_{gender_COSMIC}_{COSMIC_parameters}.pdf",
        path_plot_grade_distribution="results/plots/atacseq/exploratory/grade_distribution/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{oncotree_code_DepMap}_{gender_DepMap}_{tumor_type_GDC}_{gender_GDC}_{gender_COSMIC}_{COSMIC_parameters}.pdf",
        path_plot_age_distribution="results/plots/atacseq/exploratory/age_distribution/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{oncotree_code_DepMap}_{gender_DepMap}_{tumor_type_GDC}_{gender_GDC}_{gender_COSMIC}_{COSMIC_parameters}.pdf",
        path_plot_library_size_by_study="results/plots/atacseq/exploratory/library_size_by_study/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{oncotree_code_DepMap}_{gender_DepMap}_{tumor_type_GDC}_{gender_GDC}_{gender_COSMIC}_{COSMIC_parameters}.pdf"
    input:
        path_models_DepMap="resources/DepMap/models/{oncotree_code_DepMap}_{gender_DepMap}",
        path_manifest_GDC="resources/GDC/CNAs/{tumor_type_GDC}_{gender_GDC}_manifest",
        path_samples_COSMIC="resources/COSMIC/samples/{gender_COSMIC}_{COSMIC_parameters}",
        path_counts_atacseq="resources/atacseq/files/wCNAs/counts/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}",
        path_library_sizes="resources/atacseq/files/wCNAs/library_sizes/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}",
        path_metadata_atacseq=selected_samples_location_atacseq
    log:
        "logs/explore_atacseq_cohort/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{oncotree_code_DepMap}_{gender_DepMap}_{tumor_type_GDC}_{gender_GDC}_{gender_COSMIC}_{COSMIC_parameters}.log"
    conda:
        "../../envs/r_plotting.yml"
    shell:
        """
        Rscript workflow/scripts/explore_cohort_atacseq.R \
        --input_path_models_DepMap={input.path_models_DepMap} \
        --input_path_manifest_GDC={input.path_manifest_GDC} \
        --input_path_samples_COSMIC={input.path_samples_COSMIC} \
        --input_path_counts_atacseq={input.path_counts_atacseq} \
        --input_path_library_sizes={input.path_library_sizes} \
        --input_path_metadata_atacseq={input.path_metadata_atacseq} \
        --output_path_samples_barplot={output.path_samples_barplot} \
        --output_path_plot_samples_per_study_by_condition={output.path_plot_samples_per_study_by_condition} \
        --output_path_plot_subtype_distribution={output.path_plot_subtype_distribution} \
        --output_path_plot_grade_distribution={output.path_plot_grade_distribution} \
        --output_path_plot_age_distribution={output.path_plot_age_distribution} \
        --output_path_plot_library_size_by_study={output.path_plot_library_size_by_study} 2>{log}
        """