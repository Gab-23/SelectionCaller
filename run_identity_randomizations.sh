#!/usr/bin/env bash
#SBATCH --job-name=id_rand
#SBATCH --partition=medium
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs_pilot/identity/rand_%A_%a.out
# =====================================================================================
# Predictor identity: label randomizations = the null of the concordance test (Gab, 20/11/2026).
# One SLURM array task = one randomization: pilot_predictor_identity.py with --label_shuffle_seed = the array index,
# all requested omics in the same process (one shuffle shared by all omics). Output: $OUT/shuffle_<index>/.
# A task whose output is already complete (swap_deltas_summary.tsv present) is skipped, so a partial array can be
# resubmitted as is. Gather the results with collect_identity_randomizations.py.
#
# Launch from the root of the repository (cluster), e.g. 100 randomizations, at most 10 running at a time:
#   sbatch --array=1-100%10 --export=ALL,DS=data/feature_matrix/<dataset>,OUT=analyses/v4/BRCA/1000000bp/predictor_identity_na0_rand,BLOCK_BINS=2 run_identity_randomizations.sh
# At 20 kb: BLOCK_BINS=100 and more CPUs / memory, e.g. sbatch --cpus-per-task=16 --mem=48G --export=ALL,...,N_JOBS=16
# Overridable variables: DS (required), OUT (required), BLOCK_BINS (required: 2 at 1 Mb, 100 at 20 kb),
#   OMICS, SEEDS, N_JOBS, LABELS, NA0 (1 = --gene_na_to_zero, as the truth runs of 20/11/2026),
#   N_SHIFT (default 0: S-block shifts of the O+S model; NOTE (23/11/2026): the O+S model uses every cluster, so it does
#   not depend on the labels and the truth run already has it and its shift null; in a shuffle the shifts would only
#   cost time. The unshifted O+S fit still runs once per omic)
# The truth run must use the same DS, LABELS, SEEDS, BLOCK_BINS and NA0, otherwise truth and null are not comparable.
# =====================================================================================
# environment (before "set -u": the conda activation scripts use undefined variables)
source ~/miniforge3/etc/profile.d/conda.sh
conda activate .snakemake/conda/56a4e7b4dc24521ebcee3d03f716b893_

set -euo pipefail

: "${DS:?set DS (feature matrix)}"
: "${OUT:?set OUT (output folder of the randomizations)}"
: "${BLOCK_BINS:?set BLOCK_BINS (2 at 1 Mb, 100 at 20 kb, as the truth run)}"
: "${SLURM_ARRAY_TASK_ID:?run as an array: sbatch --array=1-100%10 ...}"
OMICS=${OMICS:-"cna_ampl cna_del atac mut_syn mut_mis mut_trunc"}
SEEDS=${SEEDS:-"1 2 3"}
N_JOBS=${N_JOBS:-${SLURM_CPUS_PER_TASK:-4}}
LABELS=${LABELS:-config/predictor_identity_labels.tsv}
NA0=${NA0:-1}
N_SHIFT=${N_SHIFT:-0}

SHUFFLE=$SLURM_ARRAY_TASK_ID
DIR=$OUT/shuffle_$SHUFFLE
if [ -s "$DIR/swap_deltas_summary.tsv" ]; then
    echo "shuffle $SHUFFLE: $DIR/swap_deltas_summary.tsv exists, skipped"
    exit 0
fi
NA0_FLAG=""
[ "$NA0" = "1" ] && NA0_FLAG="--gene_na_to_zero"

echo "commit $(git rev-parse HEAD); shuffle $SHUFFLE; dataset $DS; omics $OMICS; seeds $SEEDS; block_bins $BLOCK_BINS; NA0 $NA0; n_shift $N_SHIFT"
python -u pilot_predictor_identity.py --dataset "$DS" --labels "$LABELS" --out_dir "$DIR" \
    --omics $OMICS --seeds $SEEDS --n_jobs "$N_JOBS" --block_bins "$BLOCK_BINS" --max_buffer_bins "$BLOCK_BINS" \
    --label_shuffle_seed "$SHUFFLE" --n_shift "$N_SHIFT" $NA0_FLAG
