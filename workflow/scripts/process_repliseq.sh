
BINFILE=${20}; OUTFILE_COVERAGE=${21}

source ~/miniconda3/etc/profile.d/conda.sh # initialize conda
conda activate atacseq_preproc # activate environment
echo "Activated ATAC-seq preprocessing environment"

if [ ! -f "$NO_BLACKLISTED_SORTED_BAM.bai" ]; then
	echo "Indexing canon BAM file"
	# index BAM
	samtools index -@ $NUM_THREADS $NO_BLACKLISTED_SORTED_BAM
else
	echo "Skipping noblack BAM file indexing"
fi

if [ ! -f $OUTFILE_COVERAGE ]; then
	# compute coverage
	bedtools coverage -counts -a $BINFILE -b $NO_BLACKLISTED_SORTED_BAM > $OUTFILE_COVERAGE
else
	echo "Skipping coverage computation"
fi

echo "Done"
exit 0