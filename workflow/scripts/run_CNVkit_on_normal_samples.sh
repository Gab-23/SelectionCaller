NORMAL_SAMPLES="SRR27726627 
                SRR27726629 
                SRR27726631 
                SRR27726633 
                SRR27726636 
                SRR27726638 
                SRR27726640 
                SRR27726642 
                SRR27726644 
                SRR24401460 
                SRR24401461 
                SRR24401463 
                SRR24401465 
                SRR24401467 
                SRR24401469"
                
REGEXP=$(echo $NORMAL_SAMPLES | sed 's/ /|/g')
NORMAL_PATHS=$(ls resources/atacseq/files/bam/*/*.fixed.bam | grep -E "$REGEXP")

for NORMAL_SAMPLE in $NORMAL_SAMPLES
do
NORMAL_PATH=$(ls resources/atacseq/files/bam/*/*.fixed.bam | grep $NORMAL_SAMPLE)
ACTUAL_NORMAL_PATHS=$(echo $NORMAL_PATHS | sed 's/ /\n/g' | grep -v $NORMAL_PATH)
echo $NORMAL_SAMPLE
cnvkit.py batch $NORMAL_PATH -n $ACTUAL_NORMAL_PATHS \
        --fasta "resources/reference_data/genome/hg38.fa" \
        --method wgs \
        --target-avg-size 20000 \
        --targets "resources/backbone/hg38_binned_20000_0.filtered.sorted.bed" \
        --processes 10 \
        -d "resources/atacseq/files/cnvkit_N/"
done

TO_PROCESS=$(ls resources/atacseq/files/cnvkit_N/*_hg38_20.fixed.cns)

for FILE in $TO_PROCESS;
do
BASENAME=$(basename $FILE | cut -d "." -f 1)
cat $FILE | awk -v OFS='\t' '{{print $1,$2,$3,$5,2^$5,$6,$7,$8}}' | tail -n+2 | awk '{{print $1"\t"$2"\t"$3"\t"$5}}' > "$BASENAME".processed
done

TO_PROCESS=$(ls resources/atacseq/files/cnvkit_N/*_hg38_20.processed)

for FILE in $TO_PROCESS;
do
BASENAME=$(basename $FILE | cut -d "." -f 1)
bedtools closest -a "resources/backbone/hg38_binned_20000_0.filtered.sorted.bed" \
				 -b <(bedtools sort -i $FILE) \
				 -wa \
				 -wb | \
				 awk -v OFS='\t' '{{print $1,$2,$3,$4,$5"_"$6"_"$7,$8}}' > "$BASENAME".scaling
done
