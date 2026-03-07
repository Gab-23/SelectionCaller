#!/bin/bash

SNAKEFILE=$1
CONFIGFILE=$2
CORES=$3

source ~/anaconda3/etc/profile.d/conda.sh
conda activate snakemake-env

snakemake --snakefile $SNAKEFILE --configfile $CONFIGFILE --unlock
snakemake --snakefile $SNAKEFILE --configfile $CONFIGFILE --cores $CORES --use-conda