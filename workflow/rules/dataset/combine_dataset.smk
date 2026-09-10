rule combine_dataset:
    output:
        path_dataset="data/feature_matrix/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_{smoothing_size}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}_{gender_COSMIC}_{COSMIC_parameters}_{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}_{tissue_type_Roadmap}_{sample_type}_{liftover_str}_{file_format}_{tumor_type_IntOGen}_{haplo_thr}",
        path_corr_plot="results/plots/dataset/corr/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_{smoothing_size}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}_{gender_COSMIC}_{COSMIC_parameters}_{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}_{tissue_type_Roadmap}_{sample_type}_{liftover_str}_{file_format}_{tumor_type_IntOGen}_{haplo_thr}.pdf",
        path_NA_frac_plot="results/plots/dataset/NA/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_{smoothing_size}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}_{gender_COSMIC}_{COSMIC_parameters}_{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}_{tissue_type_Roadmap}_{sample_type}_{liftover_str}_{file_format}_{tumor_type_IntOGen}_{haplo_thr}.pdf",
        path_zero_frac_plot="results/plots/dataset/zero/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_{smoothing_size}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}_{gender_COSMIC}_{COSMIC_parameters}_{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}_{tissue_type_Roadmap}_{sample_type}_{liftover_str}_{file_format}_{tumor_type_IntOGen}_{haplo_thr}.pdf",
        path_boxplots="results/plots/dataset/boxplots/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_{smoothing_size}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}_{gender_COSMIC}_{COSMIC_parameters}_{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}_{tissue_type_Roadmap}_{sample_type}_{liftover_str}_{file_format}_{tumor_type_IntOGen}_{haplo_thr}.pdf",
        path_boxplots_scaled="results/plots/dataset/boxplots_scaled/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_{smoothing_size}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}_{gender_COSMIC}_{COSMIC_parameters}_{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}_{tissue_type_Roadmap}_{sample_type}_{liftover_str}_{file_format}_{tumor_type_IntOGen}_{haplo_thr}.pdf",
        path_boxplots_0s_non0s_muts="results/plots/dataset/boxplots_0s_non0s_muts/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_{smoothing_size}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}_{gender_COSMIC}_{COSMIC_parameters}_{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}_{tissue_type_Roadmap}_{sample_type}_{liftover_str}_{file_format}_{tumor_type_IntOGen}_{haplo_thr}.pdf",
        path_boxplots_0s_non0s_chromstates="results/plots/dataset/boxplots_0s_non0s_chromstates/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_{smoothing_size}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}_{gender_COSMIC}_{COSMIC_parameters}_{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}_{tissue_type_Roadmap}_{sample_type}_{liftover_str}_{file_format}_{tumor_type_IntOGen}_{haplo_thr}.pdf",
        path_bool_multicollinearity="results/dataset/multicollinearity/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_{smoothing_size}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}_{gender_COSMIC}_{COSMIC_parameters}_{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}_{tissue_type_Roadmap}_{sample_type}_{liftover_str}_{file_format}_{tumor_type_IntOGen}_{haplo_thr}"
    input:
        path_DARs="data/atacseq/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}",
        path_centromeres="data/distances/centromeres/{assembly_code}_{window_mode}_{window_size}_{step}",
        path_telomeres="data/distances/telomeres/{assembly_code}_{window_mode}_{window_size}_{step}",
        path_coding_mutations="data/num_overlaps/mutations/coding/{assembly_code}_{window_mode}_{window_size}_{step}_{gender_COSMIC}_{COSMIC_parameters}",
        path_noncoding_mutations="data/num_overlaps/mutations/noncoding/{assembly_code}_{window_mode}_{window_size}_{step}_{gender_COSMIC}_{COSMIC_parameters}",
        path_dist_to_essential="data/distances/essential/{assembly_code}_{window_mode}_{window_size}_{step}_{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}",
        # NOTE (fix): this used to point at resources/epigenomes/checks/processing/... - the
        # intermediate file that just contains a *path string* to the per-epigenome table (see
        # process_epigenomes.smk), not an actual binID-keyed table. merge(..., by="binID") would
        # have failed immediately (confirmed with a toy test: no "binID" column exists in that
        # file at all). This now points to aggregate_epigenomes' real, binID-keyed output.
        path_epigenomes="data/epigenomes/{tissue_type_Roadmap}_{sample_type}_{liftover_str}_{file_format}_{assembly_code}_{window_mode}_{window_size}_{step}",
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
        # NOTE (fix): added _v{gencode_version} - ohnologs/paralogs outputs were re-activated in
        # Phase 7 WITH this wildcard in their filename (their result genuinely depends on which
        # GENCODE version annotated the backbone), so this reference needed updating to match.
        path_ohnologs="data/ohnologs_and_paralogs/ohnologs/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}",
        path_paralogs="data/ohnologs_and_paralogs/paralogs/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}",
        # NOTE (fix): combine_dataset.R reads args$input_path_haploinsufficiency and merges it
        # in like every other feature, but this input was missing entirely from the rule -
        # confirmed with a toy test that the script expects it (NULL path would have crashed
        # read.delim()).
        path_haploinsufficiency="data/haploinsufficiency/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}_{haplo_thr}",
        # Added per user request: combine_dataset.R never had a Hi-C input at all - the user
        # confirmed Hi-C (and RepliSeq, before Phase 10) belongs in this feature matrix. No new
        # wildcards needed: Hi-C's own output only depends on assembly_code/mapq_thr/
        # window_mode/window_size/step, all of which this rule already has.
        path_hic="data/hic/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}.eigen"
    log:
        "logs/combine_dataset/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_{smoothing_size}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}_{gender_COSMIC}_{COSMIC_parameters}_{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}_{tissue_type_Roadmap}_{sample_type}_{liftover_str}_{file_format}_{tumor_type_IntOGen}_{haplo_thr}.log"
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
        --input_path_haploinsufficiency={input.path_haploinsufficiency} \
        --input_path_hic={input.path_hic} \
        --output_path_dataset={output.path_dataset} \
        --output_path_corr={output.path_corr_plot} \
        --output_path_NA_frac_plot={output.path_NA_frac_plot} \
        --output_path_zero_frac_plot={output.path_zero_frac_plot} \
        --output_path_boxplots={output.path_boxplots} \
        --output_path_boxplots_scaled={output.path_boxplots_scaled} \
        --output_path_boxplots_0s_non0s_muts={output.path_boxplots_0s_non0s_muts} \
        --output_path_boxplots_0s_non0s_chromstates={output.path_boxplots_0s_non0s_chromstates} \
        --output_path_bool_multicollinearity={output.path_bool_multicollinearity} 2>{log}
        """
