rule combine_dataset:
    output:
        path_dataset="data/feature_matrix/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_{smoothing_size}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}_{gender_COSMIC}_{COSMIC_parameters}_{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}_{tissue_type_Roadmap}_{sample_type}_{liftover_str}_{file_format}_{tumor_type_IntOGen}",
        path_corr_plot="results/plots/dataset/corr/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_{smoothing_size}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}_{gender_COSMIC}_{COSMIC_parameters}_{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}_{tissue_type_Roadmap}_{sample_type}_{liftover_str}_{file_format}_{tumor_type_IntOGen}",
        path_NA_plot="results/plots/dataset/NA/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_{smoothing_size}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}_{gender_COSMIC}_{COSMIC_parameters}_{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}_{tissue_type_Roadmap}_{sample_type}_{liftover_str}_{file_format}_{tumor_type_IntOGen}",
        path_boxplots="results/plots/dataset/boxplots/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_{smoothing_size}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}_{gender_COSMIC}_{COSMIC_parameters}_{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}_{tissue_type_Roadmap}_{sample_type}_{liftover_str}_{file_format}_{tumor_type_IntOGen}" 
    input:
        path_DARs="data/atacseq/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}",
        path_centromeres="data/distances/centromeres/{assembly_code}_{window_mode}_{window_size}_{step}",
        path_telomeres="data/distances/telomeres/{assembly_code}_{window_mode}_{window_size}_{step}",
        path_coding_mutations="data/num_overlaps/mutations/coding/{assembly_code}_{window_mode}_{window_size}_{step}_{gender_COSMIC}_{COSMIC_parameters}",
        path_noncoding_mutations="data/num_overlaps/mutations/noncoding/{assembly_code}_{window_mode}_{window_size}_{step}_{gender_COSMIC}_{COSMIC_parameters}",
        path_dist_to_essential="data/distances/essential/{assembly_code}_{window_mode}_{window_size}_{step}_{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}",
        path_epigenomes="resources/epigenomes/checks/processing/{tissue_type_Roadmap}_{sample_type}_{liftover_str}_{file_format}_{assembly_code}_{window_mode}_{window_size}_{step}",
        path_num_fragile_sites="data/num_overlaps/fragile_sites/{assembly_code}_{window_mode}_{window_size}_{step}",
        path_dist_fragile_sites="data/distances/fragile_sites/{assembly_code}_{window_mode}_{window_size}_{step}",
        path_GC_content="data/GC/{assembly_code}_{window_mode}_{window_size}_{step}",
        path_ampl_frequency="data/num_overlaps/CNAs/ampl/{assembly_code}_{window_mode}_{window_size}_{step}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}",
        path_del_frequency="data/num_overlaps/CNAs/del/{assembly_code}_{window_mode}_{window_size}_{step}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}",
        path_num_overlaps_protein_coding="data/num_overlaps/PC/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}",
        path_num_overlaps_lncRNA="data/num_overlaps/lncRNA/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}",
        path_num_overlaps_miRNA="data/num_overlaps/miRNA/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}",
        path_num_overlaps_rRNA="data/num_overlaps/rRNA/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}",
        path_dist_protein_coding="data/distances/PC/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}",
        path_dist_lncRNA="data/distances/lncRNA/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}",
        path_dist_miRNA="data/distances/miRNA/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}",
        path_dist_rRNA="data/distances/rRNA/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}",
        path_TSG="data/distances/TSG/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}_{tumor_type_IntOGen}",
        path_OG="data/distances/OG/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}_{tumor_type_IntOGen}",
        path_repliseq="data/repliseq/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_{smoothing_size}",
        path_ohnologs="data/ohnologs_and_paralogs/ohnologs/{assembly_code}_{window_mode}_{window_size}_{step}",
        path_paralogs="data/ohnologs_and_paralogs/paralogs/{assembly_code}_{window_mode}_{window_size}_{step}"
    log:
        "logs/combine_dataset/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_{smoothing_size}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}_{gender_COSMIC}_{COSMIC_parameters}_{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}_{tissue_type_Roadmap}_{sample_type}_{liftover_str}_{file_format}_{tumor_type_IntOGen}.log"
    conda:
        "../../envs/r_limma.yml"
    shell:
        """
        Rscript workflow/scripts/combine_dataset.R \
        --input_path_DARs={input.path_DARs} \
        --input_path_centromeres={input.path_centromeres} \
        --input_path_telomeres={input.path_telomeres} \
        --input_path_coding_mutations={input.path_coding_mutations} \
        --input_path_noncoding_mutations={input.path_noncoding_mutations} \
        --input_path_dist_to_essential={input.path_dist_to_essential} \
        --input_path_epigenomes={input.path_epigenomes} \
        --input_path_num_fragile_sites={input.path_num_fragile_sites} \
        --input_path_dist_fragile_sites={input.path_dist_fragile_sites} \
        --input_path_GC_content={input.path_GC_content} \
        --input_path_ampl_frequency={input.path_ampl_frequency} \
        --input_path_del_frequency={input.path_del_frequency} \
        --input_path_num_overlaps_protein_coding={input.path_num_overlaps_protein_coding} \
        --input_path_num_overlaps_lncRNA={input.path_num_overlaps_lncRNA} \
        --input_path_num_overlaps_miRNA={input.path_num_overlaps_miRNA} \
        --input_path_num_overlaps_rRNA={input.path_num_overlaps_rRNA} \
        --input_path_dist_protein_coding={input.path_dist_protein_coding} \
        --input_path_dist_lncRNA={input.path_dist_lncRNA} \
        --input_path_dist_miRNA={input.path_dist_miRNA} \
        --input_path_dist_rRNA={input.path_dist_rRNA} \
        --input_path_TSG={input.path_TSG} \
        --input_path_OG={input.path_OG} \
        --input_path_repliseq={input.path_repliseq} \
        --input_path_ohnologs={input.path_ohnologs} \
        --input_path_paralogs={input.path_paralogs} \
        --output_path_dataset={output.path_dataset} \
        --output_path_corr={output.path_corr_plot} \
        --output_path_NA_plot={output.path_NA_plot} \
        --output_path_boxplots={output.path_boxplots}
        """