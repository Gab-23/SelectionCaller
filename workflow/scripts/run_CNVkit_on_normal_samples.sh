#!/bin/bash
#SBATCH --job-name=CNVkit_single
#SBATCH --partition=general
#SBATCH --output=/hpcnfs/scratch/MS/ieo7429/SelectionCaller/CNVkit_%j.out
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=50
#SBATCH --mem=100G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ieo7429@ieo.it

source ~/miniforge3/etc/profile.d/conda.sh
conda activate cnvkit_env

BAM_DIR="resources/atacseq/files/bam"
DEST_DIR="resources/atacseq/files/cnvkit_N"
BACKBONE="resources/backbone/hg38_binned_20000_0.filtered.sorted.bed"
GENOME="resources/reference_data/genome/hg38.fa"
WINDOW_SIZE=20000

NORMAL_SAMPLES="SRR27726627 SRR27726629 SRR27726631 SRR27726633 
                SRR27726636 SRR27726638 SRR27726640 SRR27726642 
                SRR27726644 SRR24401460 SRR24401461 SRR24401463 
                SRR24401465 SRR24401467 SRR24401469"
                
REGEXP=$(echo $NORMAL_SAMPLES | sed 's/ /|/g')
NORMAL_PATHS=$(ls "$BAM_DIR"/*/*.fixed.bam | grep -E "$REGEXP")

for NORMAL_SAMPLE in $NORMAL_SAMPLES
do
NORMAL_PATH=$(ls "$BAM_DIR"/*/*.fixed.bam | grep $NORMAL_SAMPLE)
ACTUAL_NORMAL_PATHS=$(echo $NORMAL_PATHS | sed 's/ /\n/g' | grep -v $NORMAL_PATH)
echo $NORMAL_SAMPLE
cnvkit.py batch $NORMAL_PATH -n $ACTUAL_NORMAL_PATHS \
        --fasta $GENOME \
        --method wgs \
        --target-avg-size $WINDOW_SIZE \
        --targets $BACKBONE \
        --processes 50 \
        -d $DEST_DIR
done

TO_PROCESS=$(ls "$DEST_DIR"/*_hg38_20.fixed.cns)

for FILE in $TO_PROCESS;
do
BASENAME=$(basename $FILE | cut -d "." -f 1)
cat $FILE | \
awk -v OFS='\t' '{{print $1,$2,$3,$5,2^$5,$6,$7,$8}}' | \
tail -n+2 | \
awk '{{print $1"\t"$2"\t"$3"\t"$5}}' > \
"$DEST_DIR/$BASENAME".processed
done

conda deactivate
conda activate bedtools_env

TO_PROCESS=$(ls "$DEST_DIR"/*_hg38_20.processed)

for FILE in $TO_PROCESS;
do
BASENAME=$(basename $FILE | cut -d "." -f 1)
bedtools closest -a $BACKBONE \
		 -b <(bedtools sort -i $FILE) \
		 -wa \
		 -wb | \
		 awk -v OFS='\t' '{{print $1,$2,$3,$4,$5"_"$6"_"$7,$8}}' > \
                 "$DEST_DIR/$BASENAME".scaling
done
