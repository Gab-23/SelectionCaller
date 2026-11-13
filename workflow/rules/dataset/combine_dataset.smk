# Distance features listed in capped_distance_features (config) are read from data/distances_capped/
# (rule cap_distance_feature) instead of data/distances/. Names must match the data/distances/ subfolder.
if len(config["dist_caps"]) != 1:
    raise ValueError("dist_caps is a config-only parameter (not a wildcard): give exactly one value")

def dist_input(feature, path):
    if feature in config["capped_distance_features"]:
        return path.replace("data/distances/", "data/distances_capped/", 1)
    return path

# NOTE (29/10/2026): input del dataset ESCLUDIBILI da config (parametro solo da config, NON un wildcard):
#   dataset_excluded_inputs:
#     - DARs          # es. coorti senza ATAC tumore vs normale (polmone)
# I nomi sono quelli degli input della regola senza il prefisso "path_" (DARs, hic, repliseq, GTEx, ...). Un input escluso
# non viene richiesto (Snakemake non costruisce la sua catena) ne' passato a combine_dataset.R, e le sue colonne mancano dal
# dataset. Senza il parametro (o con lista vuota) la regola e' identica a prima.
# ATTENZIONE: non entra nel nome del dataset. Escludere un input in una coorte che ha gia' un dataset completo con lo stesso
# nome lo sovrascriverebbe con una versione senza quelle colonne: in quel caso cambiare un altro parametro del nome.
COMBINE_INPUT_NAMES = ["DARs", "centromeres", "telomeres", "coding_mutations", "coding_mutation_counts", "CDS", "coding_sites",
                       "noncoding_mutations", "dist_to_essential", "bin_level_essentiality", "epigenomes", "num_fragile_sites",
                       "dist_fragile_sites", "GC_content", "ampl_frequency", "del_frequency", "num_overlaps_protein_coding",
                       "num_overlaps_lncRNA", "num_overlaps_miRNA", "num_overlaps_rRNA", "dist_protein_coding", "dist_lncRNA",
                       "dist_miRNA", "dist_rRNA", "TSG", "OG", "repliseq", "ohnologs", "paralogs", "haploinsufficiency", "hic",
                       "interactors", "GTEx"]
DATASET_EXCLUDED_INPUTS = set(config.get("dataset_excluded_inputs", None) or [])
_unknown = DATASET_EXCLUDED_INPUTS - set(COMBINE_INPUT_NAMES)
if _unknown:
    raise ValueError(f"dataset_excluded_inputs: nomi sconosciuti {sorted(_unknown)}; validi: {COMBINE_INPUT_NAMES}")

def opt_input(name, path):
    """Percorso dell'input, oppure [] (nessun file richiesto) se l'input e' escluso da config."""
    return [] if name in DATASET_EXCLUDED_INPUTS else path

rule combine_dataset:
    output:
        path_dataset="data/feature_matrix/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_{smoothing_size}_{repliseq_cell_line}_{hic_cell_line}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}_{gender_COSMIC}_{COSMIC_parameters}_{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}_{tissue_type_Roadmap}_{sample_type}_{liftover_str}_{file_format}_{tumor_type_IntOGen}_{haplo_thr}_{GTEx_cohort}_{GTEx_release}",
        path_corr_plot="results/plots/dataset/corr/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_{smoothing_size}_{repliseq_cell_line}_{hic_cell_line}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}_{gender_COSMIC}_{COSMIC_parameters}_{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}_{tissue_type_Roadmap}_{sample_type}_{liftover_str}_{file_format}_{tumor_type_IntOGen}_{haplo_thr}_{GTEx_cohort}_{GTEx_release}.pdf",
        path_NA_frac_plot="results/plots/dataset/NA/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_{smoothing_size}_{repliseq_cell_line}_{hic_cell_line}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}_{gender_COSMIC}_{COSMIC_parameters}_{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}_{tissue_type_Roadmap}_{sample_type}_{liftover_str}_{file_format}_{tumor_type_IntOGen}_{haplo_thr}_{GTEx_cohort}_{GTEx_release}.pdf",
        path_zero_frac_plot="results/plots/dataset/zero/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_{smoothing_size}_{repliseq_cell_line}_{hic_cell_line}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}_{gender_COSMIC}_{COSMIC_parameters}_{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}_{tissue_type_Roadmap}_{sample_type}_{liftover_str}_{file_format}_{tumor_type_IntOGen}_{haplo_thr}_{GTEx_cohort}_{GTEx_release}.pdf",
        path_boxplots="results/plots/dataset/boxplots/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_{smoothing_size}_{repliseq_cell_line}_{hic_cell_line}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}_{gender_COSMIC}_{COSMIC_parameters}_{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}_{tissue_type_Roadmap}_{sample_type}_{liftover_str}_{file_format}_{tumor_type_IntOGen}_{haplo_thr}_{GTEx_cohort}_{GTEx_release}.pdf",
        path_boxplots_scaled="results/plots/dataset/boxplots_scaled/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_{smoothing_size}_{repliseq_cell_line}_{hic_cell_line}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}_{gender_COSMIC}_{COSMIC_parameters}_{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}_{tissue_type_Roadmap}_{sample_type}_{liftover_str}_{file_format}_{tumor_type_IntOGen}_{haplo_thr}_{GTEx_cohort}_{GTEx_release}.pdf",
        path_boxplots_0s_non0s_muts="results/plots/dataset/boxplots_0s_non0s_muts/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_{smoothing_size}_{repliseq_cell_line}_{hic_cell_line}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}_{gender_COSMIC}_{COSMIC_parameters}_{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}_{tissue_type_Roadmap}_{sample_type}_{liftover_str}_{file_format}_{tumor_type_IntOGen}_{haplo_thr}_{GTEx_cohort}_{GTEx_release}.pdf",
        path_boxplots_0s_non0s_chromstates="results/plots/dataset/boxplots_0s_non0s_chromstates/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_{smoothing_size}_{repliseq_cell_line}_{hic_cell_line}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}_{gender_COSMIC}_{COSMIC_parameters}_{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}_{tissue_type_Roadmap}_{sample_type}_{liftover_str}_{file_format}_{tumor_type_IntOGen}_{haplo_thr}_{GTEx_cohort}_{GTEx_release}.pdf",
        path_bool_multicollinearity="results/dataset/multicollinearity/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_{smoothing_size}_{repliseq_cell_line}_{hic_cell_line}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}_{gender_COSMIC}_{COSMIC_parameters}_{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}_{tissue_type_Roadmap}_{sample_type}_{liftover_str}_{file_format}_{tumor_type_IntOGen}_{haplo_thr}_{GTEx_cohort}_{GTEx_release}"
    input:
        path_DARs=opt_input("DARs", "data/atacseq/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}"),
        path_centromeres=opt_input("centromeres", dist_input("centromeres", "data/distances/centromeres/{assembly_code}_{window_mode}_{window_size}_{step}")),
        path_telomeres=opt_input("telomeres", dist_input("telomeres", "data/distances/telomeres/{assembly_code}_{window_mode}_{window_size}_{step}")),
        path_coding_mutations=opt_input("coding_mutations", "data/num_overlaps/mutations/coding/{assembly_code}_{window_mode}_{window_size}_{step}_{gender_COSMIC}_{COSMIC_parameters}"),
        # conteggi a livello di evento (campione x mutazione) + bp di CDS per bin (esposizione)
        path_coding_mutation_counts=opt_input("coding_mutation_counts", "data/num_overlaps/mutations/coding_counts/{assembly_code}_{window_mode}_{window_size}_{step}_{gender_COSMIC}_{COSMIC_parameters}_v{gencode_version}"),
        path_CDS=opt_input("CDS", "data/num_overlaps/CDS/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}"),
        # esposizione per classe di mutazione: sostituzioni singole possibili per bin (sites_syn / sites_mis /
        # sites_trunc / sites_other_coding), stessi trascritti dei conteggi per classe (COSMIC_transcript_mode)
        path_coding_sites=opt_input("coding_sites", "data/num_overlaps/mutations/coding_sites/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}"),
        path_noncoding_mutations=opt_input("noncoding_mutations", "data/num_overlaps/mutations/noncoding/{assembly_code}_{window_mode}_{window_size}_{step}_{gender_COSMIC}_{COSMIC_parameters}"),
        path_dist_to_essential=opt_input("dist_to_essential", dist_input("essential", "data/distances/essential/{assembly_code}_{window_mode}_{window_size}_{step}_{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}")),
        # NOTE (30/10/2026): essenzialita' DepMap a livello di bin (score medio per bin, colonne binID + ess_score), output di
        # process_bin_level_essentiality. Il path e' sotto data/num_overlaps/essential/ ma contiene SCORE, non conteggi. I bin
        # senza geni essenziali NON compaiono nel file: dopo il merge (all = TRUE) risultano NA (nessuna imputazione, da decidere).
        # Non e' una distanza: nessun dist_input / cap.
        path_bin_level_essentiality=opt_input("bin_level_essentiality", "data/num_overlaps/essential/{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}_{assembly_code}_{window_mode}_{window_size}_{step}"),
        # NOTE (fix): this used to point at resources/epigenomes/checks/processing/... - the
        # intermediate file that just contains a *path string* to the per-epigenome table (see
        # process_epigenomes.smk), not an actual binID-keyed table. merge(..., by="binID") would
        # have failed immediately (confirmed with a toy test: no "binID" column exists in that
        # file at all). This now points to aggregate_epigenomes' real, binID-keyed output.
        path_epigenomes=opt_input("epigenomes", "data/epigenomes/{tissue_type_Roadmap}_{sample_type}_{liftover_str}_{file_format}_{assembly_code}_{window_mode}_{window_size}_{step}"),
        path_num_fragile_sites=opt_input("num_fragile_sites", "data/num_overlaps/fragile_sites/{assembly_code}_{window_mode}_{window_size}_{step}"),
        path_dist_fragile_sites=opt_input("dist_fragile_sites", dist_input("fragile_sites", "data/distances/fragile_sites/{assembly_code}_{window_mode}_{window_size}_{step}")),
        path_GC_content=opt_input("GC_content", "data/GC/{assembly_code}_{window_mode}_{window_size}_{step}"),
        path_ampl_frequency=opt_input("ampl_frequency", "data/num_overlaps/CNAs/ampl/{assembly_code}_{window_mode}_{window_size}_{step}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}"),
        path_del_frequency=opt_input("del_frequency", "data/num_overlaps/CNAs/del/{assembly_code}_{window_mode}_{window_size}_{step}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}"),
        path_num_overlaps_protein_coding=opt_input("num_overlaps_protein_coding", "data/num_overlaps/PC/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}"),
        path_num_overlaps_lncRNA=opt_input("num_overlaps_lncRNA", "data/num_overlaps/lncRNA/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}"),
        path_num_overlaps_miRNA=opt_input("num_overlaps_miRNA", "data/num_overlaps/miRNA/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}"),
        path_num_overlaps_rRNA=opt_input("num_overlaps_rRNA", "data/num_overlaps/rRNA/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}"),
        path_dist_protein_coding=opt_input("dist_protein_coding", dist_input("PC", "data/distances/PC/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}")),
        path_dist_lncRNA=opt_input("dist_lncRNA", dist_input("lncRNA", "data/distances/lncRNA/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}")),
        path_dist_miRNA=opt_input("dist_miRNA", dist_input("miRNA", "data/distances/miRNA/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}")),
        path_dist_rRNA=opt_input("dist_rRNA", dist_input("rRNA", "data/distances/rRNA/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}")),
        path_TSG=opt_input("TSG", dist_input("TSG", "data/distances/TSG/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}_{tumor_type_IntOGen}")),
        path_OG=opt_input("OG", dist_input("OG", "data/distances/OG/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}_{tumor_type_IntOGen}")),
        path_repliseq=opt_input("repliseq", "data/repliseq/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_{smoothing_size}_{repliseq_cell_line}"),
        # NOTE (fix): added _v{gencode_version} - ohnologs/paralogs outputs were re-activated in
        # Phase 7 WITH this wildcard in their filename (their result genuinely depends on which
        # GENCODE version annotated the backbone), so this reference needed updating to match.
        path_ohnologs=opt_input("ohnologs", "data/ohnologs_and_paralogs/ohnologs/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}"),
        path_paralogs=opt_input("paralogs", "data/ohnologs_and_paralogs/paralogs/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}"),
        # NOTE (fix): combine_dataset.R reads args$input_path_haploinsufficiency and merges it
        # in like every other feature, but this input was missing entirely from the rule -
        # confirmed with a toy test that the script expects it (NULL path would have crashed
        # read.delim()).
        path_haploinsufficiency=opt_input("haploinsufficiency", "data/haploinsufficiency/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}_{haplo_thr}"),
        # Added per user request: combine_dataset.R never had a Hi-C input at all - the user
        # confirmed Hi-C (and RepliSeq, before Phase 10) belongs in this feature matrix. No new
        # wildcards needed: Hi-C's own output only depends on assembly_code/mapq_thr/
        # window_mode/window_size/step, all of which this rule already has.
        path_hic=opt_input("hic", "data/hic/{assembly_code}_{mapq_thr}_{window_mode}_{window_size}_{step}_{hic_cell_line}.eigen"),
        # Complex density / n. complex partners / n. PPIs in trans (CORUM + Interactome INSIDER)
        # and tissue-relative normal expression (GTEx) of the bin's genes and interactors.
        # GTEx_cohort and GTEx_release are appended to the dataset name (naming convention:
        # every parameter the matrix depends on is encoded in its name).
        path_interactors=opt_input("interactors", "data/protein_interactions/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}"),
        path_GTEx=opt_input("GTEx", "data/GTEx/features/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}_{GTEx_cohort}_{GTEx_release}")
    params:
        # NOTE (29/10/2026): argomenti --input_path_<nome>=<file> costruiti dagli input presenti (gli esclusi, vuoti, non
        # vengono passati e combine_dataset.R non li legge). Stessi nomi degli argomenti di prima.
        input_args=lambda wildcards, input: " ".join(f"--input_{name}={value}" for name, value in input.items() if value)
    log:
        "logs/combine_dataset/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{cov_thr}_{smoothing_size}_{repliseq_cell_line}_{hic_cell_line}_{minSamples_atacseq}_{minCounts_atacseq}_{minCpm_atacseq}_{SVs}_{tumor_type_GDC}_{gender_GDC}_{min_Num_Probes}_{min_abs_Segment_Mean}_{gender_COSMIC}_{COSMIC_parameters}_{oncotree_code_DepMap}_{gender_DepMap}_{essentiality_cutoff}_{num_assays_cutoff}_v{gencode_version}_{tissue_type_Roadmap}_{sample_type}_{liftover_str}_{file_format}_{tumor_type_IntOGen}_{haplo_thr}_{GTEx_cohort}_{GTEx_release}.log"
    conda:
        "../../envs/r_limma.yml"
    shell:
        """
        Rscript workflow/scripts/combine_dataset.R \
        {params.input_args} \
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