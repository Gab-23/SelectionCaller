#!/usr/bin/env bash
# Rilancia SOLO il passo 03b (driver per gene) sui run del pilota gia' completati, con la versione corretta
# di gene_level_drivers.py. Non tocca le predizioni del pilota. Uso, dalla root della repo:
#   bash rerun_gene_level.sh <cartella del run>  <dimensione del bin>
#   es. bash rerun_gene_level.sh analyses/pilot_2026-10-21_100000bp 100000
source ~/miniforge3/etc/profile.d/conda.sh
conda activate .snakemake/conda/56a4e7b4dc24521ebcee3d03f716b893_

set -euo pipefail

OUT=$1; WIN=$2; SMOOTH=$3
VER=03b_drivers_gene_level_v4

DS=data/feature_matrix/hg38_binned_${WIN}_0_20_50_${SMOOTH}_15_200_1_6_BRCA_female_50_0.2_f_breast_carcinoma_BRCA_Female_0.5_4_v45_BREAST_PrimaryCulture_hg38lift_mnemonics_BRCA_0.2_BRCA_v11
BB=resources/backbone/hg38_binned_${WIN}_0.mask.filtered.sorted.bed
ESS=resources/DepMap/essential/BRCA_Female_0.5_4
for f in "$DS" "$BB" "$ESS" gene_level_drivers.py; do [ -e "$f" ] || { echo "MANCA: $f"; exit 1; }; done

for r in chrom_full chrom_nopos genome_full genome_nopos; do
  [ -f "$OUT/02_pilot/$r/oof_predictions.tsv" ] || { echo "salto $r: predizioni assenti"; continue; }
  echo "[$(date '+%T')] $r"
  python gene_level_drivers.py --oof "$OUT/02_pilot/$r/oof_predictions.tsv" \
    --gtf resources/reference_data/gencode/gencode.v45.annotation.gtf.gz --backbone "$BB" \
    --og resources/IntOGen/OG/hg38_v45_BRCA.sorted.bed --tsg resources/IntOGen/TSG/hg38_v45_BRCA.sorted.bed \
    --cosmic resources/COSMIC/mutations/coding/f_breast_carcinoma.filtered.bed --dataset "$DS" --essential "$ESS" \
    --out_dir "$OUT/$VER/$r" > "$OUT/logs/${VER}_$r.log" 2>&1
done
{ echo; echo "== $VER ($(date '+%F %T')): gene_level_drivers.py md5 $(md5sum gene_level_drivers.py | cut -d' ' -f1)"; } >> "$OUT/MANIFEST.txt"
echo "Fatto: $OUT/$VER"