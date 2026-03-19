rule get_counts_atacseq:
    output:
        path_counts="data/atacseq/counts/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}",
        path_PCA="results/plots/atacseq/PCA/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}.pdf"
    input:
        path_check="resources/atacseq/checks/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_coverage",
        path_dir="resources/atacseq/files/coverage",
        path_metadata=atacseq_selected_samples_location
    log:
        "logs/get_counts_atacseq/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}.log"        
    conda:
        "../envs/r_backbone.yml"
    shell:
        """
        Rscript workflow/scripts/get_counts_atacseq.R \
        --input_dir={input.path_dir} \
        --input_path_metadata={input.path_metadata} \
        --minSamples_atacseq={wildcards.minSamples_atacseq} \
        --minCounts_atacseq={wildcards.minCounts_atacseq} \
        --output_path_counts={output.path_counts} \
        --output_path_PCA={output.path_PCA}
        """