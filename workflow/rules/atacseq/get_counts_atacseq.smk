rule get_counts_atacseq:
    output:
        path_library_sizes="resources/atacseq/files/wCNAs/library_sizes/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}",
        path_counts="resources/atacseq/files/wCNAs/counts/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}",
        path_model_matrix="resources/atacseq/files/wCNAs/model_matrix/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}.RData",
        path_PCA_condition="results/plots/atacseq/wCNAs/PCA/condition/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}.pdf",
        path_PCA_library_size="results/plots/atacseq/wCNAs/PCA/library_size/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}.pdf",
        path_PCA_sex="results/plots/atacseq/wCNAs/PCA/sex/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}.pdf",
        path_PCA_donor="results/plots/atacseq/wCNAs/PCA/donor/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}.pdf",
        path_PCA_study="results/plots/atacseq/wCNAs/PCA/study/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}.pdf",
        path_PCA_subtype="results/plots/atacseq/wCNAs/PCA/subtype/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}.pdf",
        path_scree_plot="results/plots/atacseq/wCNAs/scree/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}.pdf",
        path_PCA_sva_condition="results/plots/atacseq/wCNAs/PCA/condition/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_sva.pdf",
        path_PCA_sva_library_size="results/plots/atacseq/wCNAs/PCA/library_size/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_sva.pdf",
        path_PCA_sva_sex="results/plots/atacseq/wCNAs/PCA/sex/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_sva.pdf",
        path_PCA_sva_donor="results/plots/atacseq/wCNAs/PCA/donor/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_sva.pdf",
        path_PCA_sva_study="results/plots/atacseq/wCNAs/PCA/study/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_sva.pdf",
        path_PCA_sva_subtype="results/plots/atacseq/wCNAs/PCA/subtype/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_sva.pdf",
        path_scree_plot_sva="results/plots/atacseq/wCNAs/scree/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_sva.pdf",
        path_distance_plot="results/plots/atacseq/wCNAs/distance/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}.pdf",
        path_library_sizes_woCNAs="resources/atacseq/files/woCNAs/library_sizes/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}",
        path_counts_woCNAs="resources/atacseq/files/woCNAs/counts/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}",
        path_model_matrix_woCNAs="resources/atacseq/files/woCNAs/model_matrix/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}.RData",
        path_PCA_condition_woCNAs="results/plots/atacseq/woCNAs/PCA/condition/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}.pdf",
        path_PCA_library_size_woCNAs="results/plots/atacseq/woCNAs/PCA/library_size/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}.pdf",
        path_PCA_sex_woCNAs="results/plots/atacseq/woCNAs/PCA/sex/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}.pdf",
        path_PCA_donor_woCNAs="results/plots/atacseq/woCNAs/PCA/donor/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}.pdf",
        path_PCA_study_woCNAs="results/plots/atacseq/woCNAs/PCA/study/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}.pdf",
        path_PCA_subtype_woCNAs="results/plots/atacseq/woCNAs/PCA/subtype/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}.pdf",
        path_scree_plot_woCNAs="results/plots/atacseq/woCNAs/scree/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}.pdf",
        path_PCA_sva_condition_woCNAs="results/plots/atacseq/woCNAs/PCA/condition/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_sva.pdf",
        path_PCA_sva_library_size_woCNAs="results/plots/atacseq/woCNAs/PCA/library_size/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_sva.pdf",
        path_PCA_sva_sex_woCNAs="results/plots/atacseq/woCNAs/PCA/sex/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_sva.pdf",
        path_PCA_sva_donor_woCNAs="results/plots/atacseq/woCNAs/PCA/donor/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_sva.pdf",
        path_PCA_sva_study_woCNAs="results/plots/atacseq/woCNAs/PCA/study/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_sva.pdf",
        path_PCA_sva_subtype_woCNAs="results/plots/atacseq/woCNAs/PCA/subtype/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_sva.pdf",
        path_scree_plot_sva_woCNAs="results/plots/atacseq/woCNAs/scree/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_sva.pdf",
        path_distance_plot_woCNAs="results/plots/atacseq/woCNAs/distance/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}.pdf",
        path_CNVkit_lineplot="results/plots/atacseq/CNVkit/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}.pdf"
    input:
        path_check_coverage="resources/atacseq/checks/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_coverage",
        path_check_scaling_factors="resources/atacseq/checks/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_scaling_factors",
        path_metadata=selected_samples_location_atacseq,
        path_samples_to_exclude=samples_to_exclude_location_atacseq
    params:
        # bare directories, not formal Snakemake inputs (see process_CNVkit_output_atacseq.smk
        # for the same reasoning): nothing declares these whole directories as an output, since
        # they're just the accumulated per-sample files from earlier rules. The two "check"
        # inputs above already guarantee every sample's coverage/scaling-factor file exists
        # before this rule runs.
        path_dir_coverage="resources/atacseq/files/coverage",
        path_dir_scaling_factors="resources/atacseq/files/cnvkit"
    log:
        "logs/get_counts_atacseq/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}.log"        
    conda:
        "../../envs/r_limma.yml"
    shell:
        """
        Rscript workflow/scripts/get_counts_atacseq.R \
        --input_dir_coverage={params.path_dir_coverage} \
        --input_dir_scaling_factors={params.path_dir_scaling_factors} \
        --input_path_metadata={input.path_metadata} \
        --input_path_samples_to_exclude={input.path_samples_to_exclude} \
        --minSamples_atacseq={wildcards.minSamples_atacseq} \
        --minCounts_atacseq={wildcards.minCounts_atacseq} \
        --minCpm_atacseq={wildcards.minCpm_atacseq} \
        --SVs={wildcards.SVs} \
        --output_path_library_sizes={output.path_library_sizes} \
        --output_path_counts={output.path_counts} \
        --output_path_model_matrix={output.path_model_matrix} \
        --output_path_PCA_condition={output.path_PCA_condition} \
        --output_path_PCA_library_size={output.path_PCA_library_size} \
        --output_path_PCA_sex={output.path_PCA_sex} \
        --output_path_PCA_donor={output.path_PCA_donor} \
        --output_path_PCA_study={output.path_PCA_study} \
        --output_path_PCA_subtype={output.path_PCA_subtype} \
        --output_path_scree_plot={output.path_scree_plot} \
        --output_path_PCA_sva_condition={output.path_PCA_sva_condition} \
        --output_path_PCA_sva_library_size={output.path_PCA_sva_library_size} \
        --output_path_PCA_sva_sex={output.path_PCA_sva_sex} \
        --output_path_PCA_sva_donor={output.path_PCA_sva_donor} \
        --output_path_PCA_sva_study={output.path_PCA_sva_study} \
        --output_path_PCA_sva_subtype={output.path_PCA_sva_subtype} \
        --output_path_scree_plot_sva={output.path_scree_plot_sva} \
        --output_path_distance_plot={output.path_distance_plot} \
        --output_path_library_sizes_woCNAs={output.path_library_sizes_woCNAs} \
        --output_path_counts_woCNAs={output.path_counts_woCNAs} \
        --output_path_model_matrix_woCNAs={output.path_model_matrix_woCNAs} \
        --output_path_PCA_condition_woCNAs={output.path_PCA_condition_woCNAs} \
        --output_path_PCA_library_size_woCNAs={output.path_PCA_library_size_woCNAs} \
        --output_path_PCA_sex_woCNAs={output.path_PCA_sex_woCNAs} \
        --output_path_PCA_donor_woCNAs={output.path_PCA_donor_woCNAs} \
        --output_path_PCA_study_woCNAs={output.path_PCA_study_woCNAs} \
        --output_path_PCA_subtype_woCNAs={output.path_PCA_subtype_woCNAs} \
        --output_path_scree_plot_woCNAs={output.path_scree_plot_woCNAs} \
        --output_path_PCA_sva_condition_woCNAs={output.path_PCA_sva_condition_woCNAs} \
        --output_path_PCA_sva_library_size_woCNAs={output.path_PCA_sva_library_size_woCNAs} \
        --output_path_PCA_sva_sex_woCNAs={output.path_PCA_sva_sex_woCNAs} \
        --output_path_PCA_sva_donor_woCNAs={output.path_PCA_sva_donor_woCNAs} \
        --output_path_PCA_sva_study_woCNAs={output.path_PCA_sva_study_woCNAs} \
        --output_path_PCA_sva_subtype_woCNAs={output.path_PCA_sva_subtype_woCNAs} \
        --output_path_scree_plot_sva_woCNAs={output.path_scree_plot_sva_woCNAs} \
        --output_path_distance_plot_woCNAs={output.path_distance_plot_woCNAs} \
        --output_path_CNVkit_lineplot={output.path_CNVkit_lineplot} 2>{log}
        """
