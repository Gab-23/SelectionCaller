rule get_counts_atacseq:
    output:
        path_library_sizes="resources/atacseq/files/library_sizes/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}",
        path_counts="resources/atacseq/files/counts/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}",
        path_counts_sva="resources/atacseq/files/counts/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_sva",
        path_PCA_condition="results/plots/atacseq/PCA/condition/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}.pdf",
        path_PCA_library_size="results/plots/atacseq/PCA/library_size/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}.pdf",
        path_PCA_sex="results/plots/atacseq/PCA/sex/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}.pdf",
        path_PCA_donor="results/plots/atacseq/PCA/donor/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}.pdf",
        path_PCA_study="results/plots/atacseq/PCA/study/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}.pdf",
        path_PCA_subtype="results/plots/atacseq/PCA/subtype/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}.pdf",
        path_PCA_sva_condition="results/plots/atacseq/PCA/condition/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_sva.pdf",
        path_PCA_sva_library_size="results/plots/atacseq/PCA/library_size/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_sva.pdf",
        path_PCA_sva_sex="results/plots/atacseq/PCA/sex/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_sva.pdf",
        path_PCA_sva_donor="results/plots/atacseq/PCA/donor/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_sva.pdf",
        path_PCA_sva_study="results/plots/atacseq/PCA/study/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_sva.pdf",
        path_PCA_sva_subtype="results/plots/atacseq/PCA/subtype/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_sva.pdf"
    input:
        path_check="resources/atacseq/checks/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_coverage",
        path_dir="resources/atacseq/files/coverage",
        path_metadata=selected_samples_location_atacseq,
        path_samples_to_exclude=samples_to_exclude_location_atacseq
    log:
        "logs/get_counts_atacseq/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}.log"        
    conda:
        "../envs/r_limma.yml"
    shell:
        """
        Rscript workflow/scripts/get_counts_atacseq.R \
        --input_dir={input.path_dir} \
        --input_path_metadata={input.path_metadata} \
        --input_path_samples_to_exclude={input.path_samples_to_exclude} \
        --minSamples_atacseq={wildcards.minSamples_atacseq} \
        --minCounts_atacseq={wildcards.minCounts_atacseq} \
        --minCpm_atacseq={wildcards.minCpm_atacseq} \
        --output_path_library_sizes={output.path_library_sizes} \
        --output_path_counts={output.path_counts} \
        --output_path_counts_sva={output.path_counts_sva} \
        --output_path_PCA_condition={output.path_PCA_condition} \
        --output_path_PCA_library_size={output.path_PCA_library_size} \
        --output_path_PCA_sex={output.path_PCA_sex} \
        --output_path_PCA_donor={output.path_PCA_donor} \
        --output_path_PCA_study={output.path_PCA_study} \
        --output_path_PCA_subtype={output.path_PCA_subtype} \
        --output_path_PCA_sva_condition={output.path_PCA_sva_condition} \
        --output_path_PCA_sva_library_size={output.path_PCA_sva_library_size} \
        --output_path_PCA_sva_sex={output.path_PCA_sva_sex} \
        --output_path_PCA_sva_donor={output.path_PCA_sva_donor} \
        --output_path_PCA_sva_study={output.path_PCA_sva_study} \
        --output_path_PCA_sva_subtype={output.path_PCA_sva_subtype}
        """