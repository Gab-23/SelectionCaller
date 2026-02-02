rule filter_GDC_CNA_data:
    output:
        path_filtered="resources/GDC/CNAs/aggregated/{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}.filtered.bed",
        path_plot_probes="results/plots/CNAs/num_probes/{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}.pdf",
        path_plot_segment="results/plots/CNAs/segment_mean/{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}.pdf"
    input:
        "resources/GDC/CNAs/aggregated/{tumor_type_GDC}_{gender_GDC}.bed"
    log:
        "logs/filter_GDC_CNA_data/{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}.log"
    conda:
        "../envs/r_plotting.yml"
    shell:
        """
        Rscript workflow/scripts/filter_GDC_CNA_data.R \
        --input_path={input} \
        --output_path_filtered={output.path_filtered} \
        --output_path_plot_probes={output.path_plot_probes} \
        --output_path_plot_segment={output.path_plot_segment} \
        --min_Num_Probes={wildcards.min_Num_Probes} \
        --min_abs_Segment_Mean={wildcards.min_abs_Segment_Mean} 2>{log}
        """