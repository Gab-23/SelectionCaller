rule explore_cohort_atacseq:
    output:
        path_plot_samples_per_study="results/plots/atacseq/exploratory/samples_per_study/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}.pdf",
        path_plot_samples_per_condition="results/plots/atacseq/exploratory/samples_per_condition/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}.pdf",
        path_plot_subtype_distribution="results/plots/atacseq/exploratory/subtype_distribution/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}.pdf",
        path_plot_grade_distribution="results/plots/atacseq/exploratory/grade_distribution/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}.pdf",
        path_plot_samples_per_study_by_condition="results/plots/atacseq/exploratory/samples_per_study_by_condition/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}.pdf",
        path_plot_age_distribution="results/plots/atacseq/exploratory/age_distribution/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}.pdf",
        path_plot_library_size_by_study="results/plots/atacseq/exploratory/library_size_by_study/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}.pdf"
    input:
        path_metadata=selected_samples_location_atacseq,
        path_library_sizes="resources/atacseq/files/wCNAs/library_sizes/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}"
    log:
        "logs/explore_atacseq_cohort/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}.log"
    conda:
        "../../envs/r_plotting.yml"
    shell:
        """
        Rscript workflow/scripts/explore_cohort_atacseq.R \
        --input_path_metadata={input.path_metadata} \
        --input_path_library_sizes={input.path_library_sizes} \
        --output_path_plot_samples_per_study={output.path_plot_samples_per_study}\
        --output_path_plot_samples_per_condition={output.path_plot_samples_per_condition} \
        --output_path_plot_subtype_distribution={output.path_plot_subtype_distribution} \
        --output_path_plot_grade_distribution={output.path_plot_grade_distribution} \
        --output_path_plot_samples_per_study_by_condition={output.path_plot_samples_per_study_by_condition} \
        --output_path_plot_age_distribution={output.path_plot_age_distribution} \
        --output_path_plot_library_size_by_study={output.path_plot_library_size_by_study} 2>{log}
        """