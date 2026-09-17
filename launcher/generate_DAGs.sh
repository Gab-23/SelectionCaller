#!/usr/bin/env bash
# Generates one --rulegraph PDF per module of the bin-level pipeline (snakefile_dataset_bin_level).
#
# --rulegraph shows RULES and their dependencies (not concrete wildcard-resolved jobs), so it
# stays readable even for a pipeline this size - unlike --dag, which would show one node per
# job instance and, on real (non-toy) data, could mean hundreds of nodes for a single module.
#
# Each module is represented by ONE real, confirmed target path belonging to that module's
# final rule (mostly taken directly from .tests/integration - already-used, working toy
# targets; the rest checked directly against their .smk files).
#
# --rulegraph still needs to build a DAG, which means it checks that every "starting" file
# (raw/manually-downloaded data with no producing rule - fragile_sites, ohnologs,
# haploinsufficiency, DepMap, COSMIC, GDC, epigenomes) actually EXISTS on disk, even though
# it never reads its content (confirmed: MissingInputException otherwise). Rather than
# touching the real resources/ directory with placeholder files, this script runs Snakemake
# with `--directory` pointed at a throwaway scratch directory containing ONLY empty
# placeholders for those starting files - confirmed (see conversation) that `--directory`
# redirects ALL of Snakemake's relative-path resolution there, so the real repo's resources/
# is never created, read, or written by this script.
#
# Output: ONLY the final PDFs land in dags/ - no .dot, no .svg, no .log files are kept
# anywhere (snakemake's stdout is piped straight into dot -Tpdf; any stderr/error output is
# only ever printed to the console for that module, never written to a file - and the whole
# scratch directory, dot intermediates included, is removed automatically on exit via trap).
#
# Usage: bash generate_module_rulegraphs.sh [snakefile] [configfile]
# Defaults to the production Snakefile/config; pass the toy ones as args if the production
# config doesn't parse cleanly for some reason (structurally the rulegraph is the same either
# way - config values change file names, not which rules/edges appear).

set -uo pipefail  # NOT -e: one module's failure must not abort the others

SNAKEFILE="${1:-workflow/snakefile_dataset_bin_level}"
CONFIGFILE="${2:-config/config_dataset.yaml}"
OUTDIR="dags"
SCRATCH_DIR=$(mktemp -d)
trap 'rm -rf "$SCRATCH_DIR"' EXIT

mkdir -p "$OUTDIR"

if ! command -v dot >/dev/null 2>&1; then
    echo "ERROR: graphviz's 'dot' is not installed (needed to render the .pdf files)." >&2
    echo "       Install it first, e.g.: sudo apt-get install -y graphviz" >&2
    exit 1
fi

# Empty placeholders for every "starting" file (no producing rule) referenced anywhere in the
# bin-level pipeline - same set as ci.yml's dry-run-lint job, created here in $SCRATCH_DIR only,
# never in the real resources/. Content is irrelevant (--rulegraph never reads it), only
# existence and (for the two-input pairs) mtime ordering matter, and mtime ordering does NOT
# matter for --rulegraph either - it's a structural, non-executing command.
echo "Creating empty starting-file placeholders in scratch dir: $SCRATCH_DIR (real resources/ untouched)"
mkdir -p "$SCRATCH_DIR/resources/DepMap/models" "$SCRATCH_DIR/resources/DepMap/probabilities" "$SCRATCH_DIR/resources/DepMap/scores"
touch "$SCRATCH_DIR/resources/DepMap/models/Model.csv"
touch "$SCRATCH_DIR/resources/DepMap/probabilities/CRISPRGeneDependency.csv"
touch "$SCRATCH_DIR/resources/DepMap/scores/CRISPRGeneEffect.csv"
mkdir -p "$SCRATCH_DIR/resources/IntOGen/IntOGen-Drivers-20240920"
touch "$SCRATCH_DIR/resources/IntOGen/IntOGen-Drivers-20240920/Compendium_Cancer_Genes.tsv"
mkdir -p "$SCRATCH_DIR/resources/COSMIC/classification/Cosmic_Classification_Tsv_v103_GRCh38"
mkdir -p "$SCRATCH_DIR/resources/COSMIC/samples/Cosmic_Sample_Tsv_v103_GRCh38"
mkdir -p "$SCRATCH_DIR/resources/COSMIC/mutations/coding/Cosmic_GenomeScreensMutant_Tsv_v103_GRCh38"
mkdir -p "$SCRATCH_DIR/resources/COSMIC/mutations/noncoding/Cosmic_NonCodingVariants_Tsv_v103_GRCh38"
touch "$SCRATCH_DIR/resources/COSMIC/classification/Cosmic_Classification_Tsv_v103_GRCh38/Cosmic_Classification_v103_GRCh38.tsv"
touch "$SCRATCH_DIR/resources/COSMIC/samples/Cosmic_Sample_Tsv_v103_GRCh38/Cosmic_Sample_v103_GRCh38.tsv"
touch "$SCRATCH_DIR/resources/COSMIC/mutations/coding/Cosmic_GenomeScreensMutant_Tsv_v103_GRCh38/Cosmic_GenomeScreensMutant_v103_GRCh38.tsv"
touch "$SCRATCH_DIR/resources/COSMIC/mutations/noncoding/Cosmic_NonCodingVariants_Tsv_v103_GRCh38/Cosmic_NonCodingVariants_v103_GRCh38.tsv"
mkdir -p "$SCRATCH_DIR/resources/GDC/CNAs/clinical.cohort.2026-01-28" "$SCRATCH_DIR/resources/GDC/CNAs/files/dummy"
touch "$SCRATCH_DIR/resources/GDC/CNAs/clinical.cohort.2026-01-28/clinical.tsv"
touch "$SCRATCH_DIR/resources/GDC/CNAs/gdc_sample_sheet.2026-01-28.tsv"
touch "$SCRATCH_DIR/resources/GDC/CNAs/gdc_manifest.2026-01-28.170055.txt"
touch "$SCRATCH_DIR/resources/GDC/CNAs/files/dummy/dummy.seg.v2.txt"
mkdir -p "$SCRATCH_DIR/resources/fragile_sites/origin" "$SCRATCH_DIR/resources/ohnologs_and_paralogs/origin" \
         "$SCRATCH_DIR/resources/haploinsufficiency" "$SCRATCH_DIR/resources/epigenomes"
touch "$SCRATCH_DIR/resources/fragile_sites/origin/fragile_sites_li_hg38.bed"
touch "$SCRATCH_DIR/resources/fragile_sites/origin/fragile_sites_li_hg00.bed"
touch "$SCRATCH_DIR/resources/ohnologs_and_paralogs/origin/OhnologList.RData"
touch "$SCRATCH_DIR/resources/ohnologs_and_paralogs/origin/ParalogList.RData"
touch "$SCRATCH_DIR/resources/haploinsufficiency/haploinsufficiency.csv"
touch "$SCRATCH_DIR/resources/epigenomes/mapping_EID_tt_raw.tsv"

# get_counts_atacseq.smk declares two bare (non-wildcarded) DIRECTORIES as `input:` -
# resources/atacseq/files/coverage and resources/atacseq/files/cnvkit - used internally by its
# R script to glob per-sample files, but no rule in the pipeline produces exactly those bare
# paths (only wildcarded files INSIDE them, e.g. .../cnvkit/{SRA_id}/...). Confirmed by
# grepping every atacseq rule's input: block for unwildcarded paths - these are the only two.
mkdir -p "$SCRATCH_DIR/resources/atacseq/files/coverage" "$SCRATCH_DIR/resources/atacseq/files/cnvkit"

# Symlink (read-only reference, never modified) the REAL config/ directory into the scratch
# dir - config/ holds small, version-controlled files (sample lists, tissue/gender mappings,
# sequencers_info) that are read for their actual CONTENT at Snakefile parse time, before any
# rule/DAG logic runs. Confirmed directly: map_tt_acronym/map_gender/get_SRA_ids/
# get_is_patterned all read real rows, not just check existence - an empty placeholder either
# crashes outright (map_tt_acronym silently returns None on FileNotFoundError, then the 8-way
# unpack at the top of the Snakefile raises exactly the TypeError seen in this conversation)
# or, more subtly, fabricated/empty content can miss a real edge case (confirmed directly: an
# empty sequencers_info file crashes get_is_patterned with a KeyError, because pandas treats
# indexing with an empty boolean list as an empty COLUMN selector, not zero rows - a pandas
# quirk that only shows up with genuinely zero rows, which never happens with the real file).
# Reading the real config/ sidesteps all of this without maintaining a second, parallel set of
# fixture content in sync with it. This only ever READS config/ - nothing here writes to it.
ln -s "$(pwd)/config" "$SCRATCH_DIR/config"

# module_name -> one representative target path for that module's own final rule
declare -A MODULE_TARGETS=(
    [chrom_sizes]="resources/reference_data/chrom_sizes/hg00.clean.chrom.sizes"
    [blacklisted_regions]="resources/blacklist/hg00.bed"
    [backbone]="resources/backbone/hg00_binned_200_0_v1.genes_on_bins.RData"
    [gencode]="resources/reference_data/gencode/v1.sorted.bed"
    [epigenomes]="data/epigenomes/BREAST_PrimaryCulture_noliftover_mnemonics_hg00_binned_200_0"
    [IntOGen]="data/distances/TSG/hg00_binned_200_0_v1_BRCA"
    [genes]="data/distances/PC/hg00_binned_200_0_v1"
    [GC]="data/GC/hg00_binned_200_0"
    [centromeres_and_telomeres]="data/distances/centromeres/hg00_binned_200_0"
    [DepMap]="data/distances/essential/hg00_binned_200_0_BRCA_Female_0.5_1_v1"
    [gene_density]="resources/num_overlaps/density/essential/hg00_binned_200_0_BRCA_Female_0.5_1_v1_BRCA"
    [haploinsufficiency]="data/haploinsufficiency/hg00_binned_200_0_v1_0.9"
    [COSMIC]="data/num_overlaps/mutations/coding/hg00_binned_200_0_f_breast_carcinoma"
    [GDC_CNA]="data/num_overlaps/CNAs/ampl/hg00_binned_200_0_BRCA_female_50_0.2"
    [fragile_sites]="data/distances/fragile_sites/hg00_binned_200_0"
    [ohnologs_and_paralogs]="data/ohnologs_and_paralogs/ohnologs/hg00_binned_200_0_v1"
    [genome]="resources/reference_data/genome/hg00.fa.bwt"
    [repliseq]="data/repliseq/hg00_binned_200_0_20_0_500000"
    [atacseq]="data/atacseq/hg00_binned_200_0_20_2_1_1_1_BRCA_female_50_0.2"
    [HiC]="data/hic/hg00_20_binned_200_0.eigen"
    [dataset]="data/feature_matrix/hg00_binned_200_0_20_0_500000_2_1_1_1_BRCA_female_50_0.2_f_breast_carcinoma_BRCA_Female_0.5_1_v1_BREAST_PrimaryCulture_noliftover_mnemonics_BRCA_0.9"
)

FAILED=()
OK=()

for module in "${!MODULE_TARGETS[@]}"; do
    target="${MODULE_TARGETS[$module]}"
    out="$OUTDIR/${module}.rulegraph.pdf"
    dot_tmp="$SCRATCH_DIR/${module}.dot"       # scratch only - removed with $SCRATCH_DIR on exit
    err_tmp="$SCRATCH_DIR/${module}.err"        # scratch only - never kept, printed then discarded

    echo "=== $module -> $target ==="
    if snakemake -s "$(pwd)/$SNAKEFILE" --configfile "$(pwd)/$CONFIGFILE" --directory "$SCRATCH_DIR" \
        --rulegraph dot "$target" > "$dot_tmp" 2>"$err_tmp"; then
        if dot -Tpdf "$dot_tmp" > "$out" 2>"$err_tmp"; then
            echo "  [OK] $out"
            OK+=("$module")
        else
            echo "  [FAIL] dot rendering failed:"
            sed 's/^/    /' "$err_tmp"
            FAILED+=("$module")
        fi
    else
        echo "  [FAIL] snakemake --rulegraph failed:"
        sed 's/^/    /' "$err_tmp"
        FAILED+=("$module")
    fi
done

echo
echo "=== Summary ==="
echo "OK (${#OK[@]}): ${OK[*]}"
if [ "${#FAILED[@]}" -gt 0 ]; then
    echo "FAILED (${#FAILED[@]}): ${FAILED[*]}"
    echo "(error details were printed above for each failed module - nothing was written to disk for them)"
    exit 1
fi
echo "All module rulegraphs generated as PDF in $OUTDIR/ - real resources/ was never touched (used $SCRATCH_DIR)"