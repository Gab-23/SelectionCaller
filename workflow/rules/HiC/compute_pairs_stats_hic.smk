rule compute_pairs_stats_hic:
	output:
		"resources/hic/files/pairs/{SRA_id}/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}.filtered.dedup.clipped.stats"
	input:
		"resources/hic/files/pairs/{SRA_id}/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}.filtered.dedup.sorted.clipped.pairs"
	log:
		"logs/compute_pairs_stat_hic/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{SRA_id}.log"
	conda:
		"../../envs/hic_processing_env.yml"
	shell:
		"""
		pairtools stats {input} -o {output} 2>{log}
		"""

