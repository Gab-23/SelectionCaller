#!/usr/bin/env bash
#SBATCH --job-name=pilot_analysis
#SBATCH --partition=medium
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=pilot_analysis_%j.out
# =====================================================================================
# Analisi completa del pilota, da zero, in una cartella nuova e datata (21/10/2026).
#   01 tetti di rumore (CNA ampl/del dai file GDC per campione; mutazioni con e senza driver)
#   02 pilota: 2 run GENOME-WIDE (con / senza distanze da telomero e centromero); i run per cromosoma sono
#      stati tolti (22/10/2026): a 500 kb e 1 Mb i cromosomi hanno troppi pochi bin e venivano saltati in silenzio,
#      tutti con --loss glm (CNA binomiale, mutazioni Poisson con offset, ATAC gaussiana)
#   03 controllo dei driver per ciascun run
#   04 figure riassuntive
# Tutto e' registrato in MANIFEST.txt (dataset e script con md5, commit git, parametri, ambiente).
#   03b controllo dei driver PER GENE (ogni gene conta una volta) per ciascun run
#   Il pilota calcola anche il PAVIMENTO (N_FLOOR ripetizioni con il target spostato circolarmente).
# Lanciare dalla root della repo, con l'ambiente py_ml_env attivo. Una risoluzione per job:
#   sbatch --export=ALL,WIN=20000 run_pilot_analysis.sh
#   sbatch --export=ALL,WIN=100000 run_pilot_analysis.sh      (ecc.)
#   sbatch --export=ALL,WIN=1000000,SMOOTH=4000000 run_pilot_analysis.sh
# Variabili sovrascrivibili: WIN, SMOOTH, DS, BB, OUT, N_JOBS, N_FLOOR, FLOOR_SCOPE
# =====================================================================================
# ambiente (prima di "set -u": gli script di attivazione di conda usano variabili non definite)
source ~/miniforge3/etc/profile.d/conda.sh
conda activate .snakemake/conda/56a4e7b4dc24521ebcee3d03f716b893_

set -euo pipefail

WIN=${WIN:-20000}                                   # dimensione del bin (bp)
SMOOTH=${SMOOTH:-500000}                            # smoothing del dataset: 500000 per 20 kb / 100 kb, 4000000 per 500 kb / 1 Mb
DS=${DS:-data/feature_matrix/hg38_binned_${WIN}_0_20_50_${SMOOTH}_15_200_1_6_BRCA_female_50_0.2_f_breast_carcinoma_BRCA_Female_0.5_4_v45_BREAST_PrimaryCulture_hg38lift_mnemonics_BRCA_0.2_BRCA_v11}
BB=${BB:-resources/backbone/hg38_binned_${WIN}_0.mask.filtered.sorted.bed}
GTF=${GTF:-resources/reference_data/gencode/gencode.v45.annotation.gtf.gz}
OG=${OG:-resources/IntOGen/OG/hg38_v45_BRCA.sorted.bed}
TSG=${TSG:-resources/IntOGen/TSG/hg38_v45_BRCA.sorted.bed}
ESS=${ESS:-resources/DepMap/essential/BRCA_Female_0.5_4}
# blocchi della CV e tetto del cuscinetto: stessa scala FISICA a ogni risoluzione (2 Mb), minimo 2 bin
BLOCK_BINS=$(( 2000000 / WIN > 2 ? 2000000 / WIN : 2 ))
MAX_BUF_BINS=$BLOCK_BINS
N_FLOOR=${N_FLOOR:-3}
FLOOR_SCOPE=${FLOOR_SCOPE:-chrom}
GDC_FILES=${GDC_FILES:-"resources/GDC/CNAs/files/*/*.seg.v2.txt"}
COSMIC=${COSMIC:-resources/COSMIC/mutations/coding/f_breast_carcinoma.filtered.bed}
GDC_AMPL=resources/GDC/CNAs/ampl/BRCA_female_50_0.2.bed
GDC_DEL=resources/GDC/CNAs/del/BRCA_female_50_0.2.bed
OUT=${OUT:-analyses/pilot_$(date +%F)_${WIN}bp}
N_JOBS=${N_JOBS:-${SLURM_CPUS_PER_TASK:-$(nproc)}}
SEED=1
PY=${PY:-python}
SCRIPTS="pilot_omics_background.py noise_ceiling.py residual_vs_raw_drivers.py gene_level_drivers.py summary_figures_pilot.py"

log() { echo "[$(date '+%F %T')] $*"; }
for f in $DS $BB $COSMIC $GDC_AMPL $GDC_DEL $GTF $OG $TSG $ESS workflow/scripts/helpers/INTERLEAVED_KFOLD_IDXS.py $SCRIPTS; do [ -e "$f" ] || { echo "MANCA: $f"; exit 1; }; done
if [ -e "$OUT" ]; then echo "La cartella $OUT esiste gia': scegli un altro OUT (non sovrascrivo)."; exit 1; fi
mkdir -p "$OUT"/{01_noise_ceiling,02_pilot,03_drivers,03b_drivers_gene_level,04_figures,logs}

# colonne di posizione da escludere nei run "nopos" (lette dal dataset, non scritte a mano)
POS=$(head -1 "$DS" | tr '\t' '\n' | grep -iE "telomer|centromer" | tr '\n' ' ')
[ -n "$POS" ] || { echo "Nessuna colonna di telomero/centromero trovata nel dataset"; exit 1; }

# ---------------- MANIFEST ----------------
{
  echo "Analisi pilota - $(date '+%F %T')  host=$(hostname)  job=${SLURM_JOB_ID:-interattivo}"
  echo "Cartella: $OUT"
  echo "Commit git: $(git rev-parse --verify -q HEAD || echo n/a)  modifiche non committate: $(git status --porcelain 2>/dev/null | wc -l)"
  echo; echo "== Input (md5) =="
  md5sum "$DS" "$BB" "$COSMIC" "$GDC_AMPL" "$GDC_DEL" "$GTF" "$OG" "$TSG" "$ESS"
  echo "File GDC per campione: $(ls $GDC_FILES | wc -l)"
  echo; echo "== Script (md5) =="; md5sum $SCRIPTS
  echo; echo "== Parametri =="
  echo "bin=$WIN bp  smoothing=$SMOOTH bp  loss=glm  seed=$SEED  n_jobs=$N_JOBS  block_bins=$BLOCK_BINS  max_buffer_bins=$MAX_BUF_BINS"
  echo "pavimento: $N_FLOOR ripetizioni, shift circolare scope=$FLOOR_SCOPE"
  echo "colonne di posizione escluse nei run nopos: $POS"
  echo; echo "== Ambiente =="; $PY --version 2>&1
  $PY -c "import numpy, pandas, xgboost, sklearn, matplotlib; print('numpy', numpy.__version__, '| pandas', pandas.__version__, '| xgboost', xgboost.__version__, '| sklearn', sklearn.__version__, '| matplotlib', matplotlib.__version__)"
} > "$OUT/MANIFEST.txt"
cp $SCRIPTS "$OUT/logs/"          # copia esatta degli script usati
log "Manifest scritto in $OUT/MANIFEST.txt"

# ---------------- 01 tetti di rumore ----------------
for d in ampl del; do
  log "01 tetto CNA $d"
  $PY noise_ceiling.py --kind cna --events "$GDC_FILES" --direction $d --backbone "$BB" --dataset "$DS" \
      --seed $SEED --out_dir "$OUT/01_noise_ceiling/cna_$d" > "$OUT/logs/01_ceiling_cna_$d.log" 2>&1
done
log "01 tetto mutazioni (con e senza driver)"
$PY noise_ceiling.py --kind mut --events "$COSMIC" --backbone "$BB" --dataset "$DS" --exclude_drivers \
    --seed $SEED --out_dir "$OUT/01_noise_ceiling/mut" > "$OUT/logs/01_ceiling_mut.log" 2>&1

# controllo di coerenza: i segmenti selezionati devono coincidere con i file ampl/del della pipeline
{
  echo; echo "== Controlli di coerenza dei tetti =="
  for d in ampl del; do
    n_ceil=$(grep -o "segmenti dopo i filtri ($d): [0-9]*" "$OUT/logs/01_ceiling_cna_$d.log" | grep -o "[0-9]*$")
    f=$([ $d = ampl ] && echo $GDC_AMPL || echo $GDC_DEL); n_pipe=$(wc -l < "$f")
    echo "CNA $d: segmenti nel tetto = $n_ceil; righe nel file della pipeline = $n_pipe  -> $([ "$n_ceil" = "$n_pipe" ] && echo OK || echo DIVERSI)"
    grep "campioni) letti" "$OUT/logs/01_ceiling_cna_$d.log" | sed 's/^/  /'
  done
} >> "$OUT/MANIFEST.txt"

# ---------------- 02 pilota: 4 run ----------------
run_pilot() {   # nome, mode, extra
  local name=$1 mode=$2; shift 2
  log "02 pilota $name"
  $PY pilot_omics_background.py --dataset "$DS" --loss glm --mode $mode --n_jobs $N_JOBS --seed $SEED \
      --bin_size $WIN --block_bins $BLOCK_BINS --max_buffer_bins $MAX_BUF_BINS \
      --n_floor $N_FLOOR --floor_scope $FLOOR_SCOPE \
      --out_dir "$OUT/02_pilot/$name" "$@" > "$OUT/logs/02_pilot_$name.log" 2>&1
}
RUNS="genome_full genome_nopos"
run_pilot genome_full  genome
run_pilot genome_nopos genome --exclude_features $POS

# ---------------- 03 driver ----------------
for r in $RUNS; do
  log "03 driver $r"
  $PY residual_vs_raw_drivers.py --oof "$OUT/02_pilot/$r/oof_predictions.tsv" --dataset "$DS" --seed $SEED \
      --out_dir "$OUT/03_drivers/$r" > "$OUT/logs/03_drivers_$r.log" 2>&1
done

# ---------------- 03b driver per gene ----------------
for r in $RUNS; do
  log "03b driver per gene $r"
  $PY gene_level_drivers.py --oof "$OUT/02_pilot/$r/oof_predictions.tsv" --gtf "$GTF" --backbone "$BB" \
      --og "$OG" --tsg "$TSG" --cosmic "$COSMIC" --dataset "$DS" --essential "$ESS" --seed $SEED \
      --out_dir "$OUT/03b_drivers_gene_level/$r" > "$OUT/logs/03b_gene_level_$r.log" 2>&1
done

# ---------------- 04 figure ----------------
log "04 figure"
$PY summary_figures_pilot.py --dataset "$DS" \
  --run "$OUT/02_pilot/genome_full"  "Genome-wide, with chr. coordinates" \
  --run "$OUT/02_pilot/genome_nopos" "Genome-wide, no chr. coordinates" \
  --modules_run "$OUT/02_pilot/genome_nopos" \
  --ceiling cna_ampl "$OUT/01_noise_ceiling/cna_ampl/noise_ceiling.tsv" \
  --ceiling cna_del  "$OUT/01_noise_ceiling/cna_del/noise_ceiling.tsv" \
  --ceiling mut      "$OUT/01_noise_ceiling/mut/noise_ceiling.tsv" \
  --drivers_tsv "$OUT/03_drivers/genome_nopos/residual_vs_raw.tsv" \
  --drivers_run "$OUT/02_pilot/genome_nopos" \
  --out_dir "$OUT/04_figures" > "$OUT/logs/04_figures.log" 2>&1

# riepilogo finale nel manifest: famiglie usate e scala delle mutazioni riconosciuta
{
  echo; echo "== Famiglie per omica e run (dalle predizioni OOF) =="
  for r in $RUNS; do
    echo "$r: $(cut -f1,8 "$OUT/02_pilot/$r/oof_predictions.tsv" | tail -n+2 | sort -u | tr '\t' '=' | tr '\n' ' ')"
  done
  echo; echo "== Pavimento (dal log del pilota) =="
  for r in $RUNS; do
    echo "$r:"; [ -f "$OUT/02_pilot/$r/floor.tsv" ] && $PY -c "import pandas as pd; F=pd.read_csv('$OUT/02_pilot/$r/floor.tsv',sep='\t'); print(F.groupby('omic')[['R2','dev_expl']].mean().round(3).to_string())" | sed 's/^/  /' || echo "  ASSENTE"
  done
  echo; echo "== Log delle figure =="; grep -i "mutazioni: target\|pavimento" "$OUT/logs/04_figures.log" || true
  echo; echo "Completato: $(date '+%F %T')"
} >> "$OUT/MANIFEST.txt"
log "Fatto. Risultati in $OUT (vedi MANIFEST.txt)"