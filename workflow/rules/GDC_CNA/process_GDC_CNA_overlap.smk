rule process_GDC_CNA_overlap:
    output:
        path_ampl="data/num_overlaps/CNAs/ampl/{assembly_code}_{window_mode}_{window_size}_{step}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}",
        path_del="data/num_overlaps/CNAs/del/{assembly_code}_{window_mode}_{window_size}_{step}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}",
        path_plot_ampl="results/plots/CNAs/landscape/ampl/{assembly_code}_{window_mode}_{window_size}_{step}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}.pdf",
        path_plot_del="results/plots/CNAs/landscape/del/{assembly_code}_{window_mode}_{window_size}_{step}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}.pdf"
    input:
        path_ampl="resources/num_overlaps/CNAs/ampl/{assembly_code}_{window_mode}_{window_size}_{step}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}",
        path_del="resources/num_overlaps/CNAs/del/{assembly_code}_{window_mode}_{window_size}_{step}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}",
        path_manifest="resources/GDC/CNAs/{tumor_type_GDC}_{gender_GDC}_manifest"
    log:
        "logs/process_GDC_CNA_overlap/{assembly_code}_{window_mode}_{window_size}_{step}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}.log"
    conda:
        "../../envs/r_plotting.yml"
    shell:
        """
        NUM_PATIENTS=$(tail -n+2 {input.path_manifest} | wc -l)
        Rscript workflow/scripts/process_GDC_CNA_overlap.R \
        --input_path_ampl={input.path_ampl} \
        --input_path_del={input.path_del} \
        --output_path_ampl={output.path_ampl} \
        --output_path_del={output.path_del} \
        --output_path_plot_ampl={output.path_plot_ampl} \
        --output_path_plot_del={output.path_plot_del} \
        --num_patients=$NUM_PATIENTS 2> {log}
        """
