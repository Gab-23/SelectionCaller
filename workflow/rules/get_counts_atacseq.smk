rule get_counts_atacseq:
    output:
        path_counts="data/atacseq/counts/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atac}_{minCounts_atac}",
        path_PCA="results/plots/atacseq/PCA/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atac}_{minCounts_atac}.pdf"
    input:
        path_check="resources/atacseq/checks/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_coverage",
        path_dir="resources/atacseq/files/coverage"
    log:
        "logs/get_counts_atacseq/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}.log"        
    conda:
        "../envs/r_backbone.yml"
    shell:
        """
        Rscript workflow/scripts/get_counts_atacseq.R \
        --input_dir={input.path_dir} \
        --minSamples_atac={wildcards.minSamples_atac} \
        --minCounts_atac={wildcards.minCounts_atac} \
        --output_path_counts={output.path_counts} \
        --output_path_PCA={output.path_PCA}
        """