rule compute_pairs_stats_hic:
	output:
		"resources/hic/files/pairs/{SRA_id}/{assembly_code}_{mapq_thr}.filtered.dedup.stats"
	input:
		"resources/hic/files/pairs/{SRA_id}/{assembly_code}_{mapq_thr}.filtered.dedup.sorted.pairs"
	log:
		"logs/compute_pairs_stat_hic/{assembly_code}_{mapq_thr}_{SRA_id}.log"
	conda:
		"../envs/hic_processing_env.yaml"
	shell:
		"""
		pairtools stats {input} -o {output} 2>{log}
		"""

