rule deduplicate_pairs_hic:
	output:
		pairs="resources/hic/files/pairs/{SRA_id}/{assembly_code}.dedup.sorted.pairs",
		stats="resources/hic/files/pairs/{SRA_id}/{assembly_code}.dedup.stats"
	input:
		 "resources/hic/files/pairs/{SRA_id}/{assembly_code}.sorted.pairs"
	log:
		"logs/deduplicate_pairs_hic/{assembly_code}_{SRA_id}.log"
	params:
		num_cores_hic=num_cores_hic
	conda:
		"../../envs/hic_processing_env.yaml"
	shell:
		"""
		pairtools dedup \
		--nproc-in {params.num_cores_hic} \
		--nproc-out {params.num_cores_hic} \
		--max-mismatch 3 \
		--mark-dups \
		--output {output.pairs} \
		--output-stats {output.stats} \
		{input} &

		DEDUP_PID=$!
		while kill -0 $DEDUP_PID ; do
			echo "Process is still active..."
		sleep 1
		done
		"""