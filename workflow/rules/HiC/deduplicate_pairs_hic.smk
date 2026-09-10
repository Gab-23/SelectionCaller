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
		"../../envs/hic_processing_env.yml"
	shell:
		# NOTE (fix): kept the original background+heartbeat pattern (confirmed necessary by
		# the user for their environment), but added the exit-status check that was missing -
		# the original silently reported success even if `pairtools dedup` crashed in the
		# background, since neither `set -e` nor the polling loop ever inspected its exit code
		# (confirmed with a direct test). `wait $DEDUP_PID` now captures the real exit code and
		# fails the rule if it's non-zero.
		"""
		pairtools dedup \
		--nproc-in {params.num_cores_hic} \
		--nproc-out {params.num_cores_hic} \
		--max-mismatch 3 \
		--mark-dups \
		--output {output.pairs} \
		--output-stats {output.stats} \
		{input} 2>{log} &

		DEDUP_PID=$!
		while kill -0 $DEDUP_PID 2>/dev/null ; do
			echo "Process is still active..."
		sleep 1
		done

		wait $DEDUP_PID
		"""
