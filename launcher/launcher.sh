#!/bin/bash
#SBATCH --job-name=snakemake_single
#SBATCH --partition=general
#SBATCH --output=/hpcnfs/scratch/MS/ieo7429/SelectionCaller/nextflow_%j.out
#SBATCH --error=/hpcnfs/scratch/MS/ieo7429/SelectionCaller/nextflow_%j.err
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=500G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ieo7429@ieo.it

SNAKEFILE=$1
CONFIGFILE=$2
CORES=$3

SOURCE="ieo7429@146.247.19.84"
LOCAL="ieo7429@login.ieo.it"

SOURCE_BASE="/home/ieo7429/Desktop/SelectionCaller"
SOURCE_RESOURCES="$SOURCE_BASE/resources"
SOURCE_DATA="$SOURCE_BASE/data"
SOURCE_RESULTS="$SOURCE_BASE/results"
SOURCE_LOGS="$SOURCE_BASE/logs"

LOCAL_BASE="/hpcnfs/scratch/MS/ieo7429/SelectionCaller"
LOCAL_RESOURCES="$LOCAL_BASE/resources"
LOCAL_DATA="$LOCAL_BASE/data"
LOCAL_RESULTS="$LOCAL_BASE/results"
LOCAL_LOGS="$LOCAL_BASE/logs"

echo "=================================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo "=================================================="

echo "[$(date)] Pulling input data from $SOURCE..."

rsync -avz --progress \
  "$SOURCE:$SOURCE_RESOURCES" \
  "$LOCAL_BASE/" \
  || { echo "ERROR: rsync of resources failed"; exit 1; }

echo "[$(date)] Data transfer complete."

source ~/miniforge3/etc/profile.d/conda.sh
conda activate snakemake-env

snakemake --snakefile $SNAKEFILE --configfile $CONFIGFILE --unlock
snakemake --snakefile $SNAKEFILE --configfile $CONFIGFILE --cores $CORES --use-conda

EXIT_STATUS=$?

if [ $EXIT_STATUS -eq 0 ]; then
  echo "[$(date)] Pipeline finished successfully. Sending results to $SOURCE..."

  ssh "$SOURCE" "mkdir -p $SOURCE_DATA $SOURCE_RESULTS $SOURCE_LOGS"

  rsync -avz --progress \
    "$LOCAL_DATA/" \
    "$SOURCE:$SOURCE_DATA" \
    && echo "[$(date)] Data/ transferred successfully." \
    || { echo "ERROR: rsync of Data/ failed. NOT deleting local copy."; exit 1; }

  rsync -avz --progress \
    "$LOCAL_RESULTS/" \
    "$SOURCE:$SOURCE_RESULTS" \
    && echo "[$(date)] Results/ transferred successfully." \
    || { echo "ERROR: rsync of Results/ failed. NOT deleting local copy."; exit 1; }

  rsync -avz --progress \
    "$LOCAL_LOGS/" \
    "$SOURCE:$SOURCE_LOGS" \
    && echo "[$(date)] Logs/ transferred successfully." \
    || { echo "ERROR: rsync of Logs/ failed. NOT deleting local copy."; exit 1; }

  #echo "[$(date)] Cleaning up local data and results..."
  #rm -rf "$LOCAL_DATA" "$LOCAL_RESULTS" "$LOCAL_LOGS"
  #echo "[$(date)] Cleanup done."

else
  echo "ERROR: Pipeline failed (exit $EXIT_STATUS). Skipping result transfer and cleanup."
fi

echo "=================================================="
echo "End time: $(date)"
echo "Exit status: $EXIT_STATUS"
echo "=================================================="

exit $EXIT_STATUS