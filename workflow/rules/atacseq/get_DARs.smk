rule get_DARs:
    output:
        path_DARs="resources/atacseq/files/DARs/wCNAs/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}",
        path_DARs_woCNAs="resources/atacseq/files/DARs/woCNAs/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}",
        path_diagnostics_plot="results/plots/atacseq/wCNAs/diagnostics/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}.pdf",
        path_diagnostics_plot_woCNAs="results/plots/atacseq/woCNAs/diagnostics/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}.pdf",
        path_volcano_plot="results/plots/atacseq/wCNAs/volcano/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}.pdf",
        path_volcano_plot_woCNAs="results/plots/atacseq/woCNAs/volcano/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}.pdf",
        path_landscape_plot="results/plots/atacseq/wCNAs/landscape/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}.pdf",
        path_landscape_plot_woCNAs="results/plots/atacseq/woCNAs/landscape/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}.pdf",
        path_CNA_comparison_plot="results/plots/atacseq/wCNAs/comparison/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}.pdf",
        path_CNA_comparison_plot_woCNAs="results/plots/atacseq/woCNAs/comparison/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}.pdf",
        # NEW: GC-bias diagnostic (logFC / logCPM by GC-quantile bin)
        path_GC_bias_plot="results/plots/atacseq/wCNAs/GC_bias/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}.pdf",
        path_GC_bias_plot_woCNAs="results/plots/atacseq/woCNAs/GC_bias/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}.pdf",
        # NEW: per-sample GC curves (colour = condition, columns = study_id)
        path_GC_curves_plot="results/plots/atacseq/wCNAs/GC_curves/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}.pdf",
        path_GC_curves_plot_woCNAs="results/plots/atacseq/woCNAs/GC_curves/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}.pdf"
    input:
        path_metadata=selected_samples_location_atacseq,
        path_counts="resources/atacseq/files/wCNAs/counts/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}",
        path_counts_woCNAs="resources/atacseq/files/woCNAs/counts/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}",
        path_library_sizes="resources/atacseq/files/wCNAs/library_sizes/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}",
        path_library_sizes_woCNAs="resources/atacseq/files/woCNAs/library_sizes/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}",
        path_model_matrix="resources/atacseq/files/wCNAs/model_matrix/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}.RData",
        path_model_matrix_woCNAs="resources/atacseq/files/woCNAs/model_matrix/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}.RData",
        path_ampl="data/num_overlaps/CNAs/ampl/{assembly_code}_{window_mode}_{window_size}_{step}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}",
        path_del="data/num_overlaps/CNAs/del/{assembly_code}_{window_mode}_{window_size}_{step}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}",
        # NEW: per-bin GC content, same binning as the counts (output of process_GC_content)
        path_GC="data/GC/{assembly_code}_{window_mode}_{window_size}_{step}"
    params:
        # NEW: GC correction of the DA analysis (see workflow/scripts/helpers/compute_GC_offset.R).
        # Defaults keep the previous behaviour if the keys are absent from the config.
        # NOTE: these are params, not wildcards -> changing them re-runs this rule and OVERWRITES
        # the DARs (and everything downstream becomes out of date).
        GC_correction=config.get("GC_correction_atacseqs", ["none"])[0],
        GC_span=config.get("GC_span_atacseqs", [0.3])[0],
        GC_n_fit=config.get("GC_n_fit_atacseqs", [20000])[0]
    log:
        "logs/get_DARs/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}.log"
    conda:
        "../../envs/r_limma.yml"
    shell:
        """
        Rscript workflow/scripts/get_DARs.R \
        --input_path_metadata={input.path_metadata} \
        --input_path_counts={input.path_counts} \
        --input_path_counts_woCNAs={input.path_counts_woCNAs} \
        --input_path_library_sizes={input.path_library_sizes} \
        --input_path_library_sizes_woCNAs={input.path_library_sizes_woCNAs} \
        --input_path_model_matrix={input.path_model_matrix} \
        --input_path_model_matrix_woCNAs={input.path_model_matrix_woCNAs} \
        --input_path_ampl={input.path_ampl} \
        --input_path_del={input.path_del} \
        --input_path_GC={input.path_GC} \
        --GC_correction={params.GC_correction} \
        --GC_span={params.GC_span} \
        --GC_n_fit={params.GC_n_fit} \
        --output_path_DARs={output.path_DARs} \
        --output_path_DARs_woCNAs={output.path_DARs_woCNAs} \
        --output_path_diagnostics_plot={output.path_diagnostics_plot} \
        --output_path_diagnostics_plot_woCNAs={output.path_diagnostics_plot_woCNAs} \
        --output_path_volcano_plot={output.path_volcano_plot} \
        --output_path_volcano_plot_woCNAs={output.path_volcano_plot_woCNAs} \
        --output_path_landscape_plot={output.path_landscape_plot} \
        --output_path_landscape_plot_woCNAs={output.path_landscape_plot_woCNAs} \
        --output_path_CNA_comparison_plot={output.path_CNA_comparison_plot} \
        --output_path_CNA_comparison_plot_woCNAs={output.path_CNA_comparison_plot_woCNAs} \
        --output_path_GC_bias_plot={output.path_GC_bias_plot} \
        --output_path_GC_bias_plot_woCNAs={output.path_GC_bias_plot_woCNAs} \
        --output_path_GC_curves_plot={output.path_GC_curves_plot} \
        --output_path_GC_curves_plot_woCNAs={output.path_GC_curves_plot_woCNAs} 2>{log}
        """