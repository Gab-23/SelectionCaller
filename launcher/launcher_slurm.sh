#!/bin/bash

#SBATCH --job-name=snakemake_single
#SBATCH --partition=general
#SBATCH --cpus-per-task=20
#SBATCH --mem=500G
#SBATCH --time=96:00:00
#SBATCH --output=snakemake_%j.out

source ~/miniforge3/etc/profile.d/conda.sh
conda activate snakemake-env

# Run locally using the allocated resources
snakemake --configfile config/config.yaml --snakefile workflow/workflow_bin_level --unlock
snakemake --configfile config/config.yaml --snakefile workflow/workflow_bin_level --use-conda --rerun-incomplete --cores 10