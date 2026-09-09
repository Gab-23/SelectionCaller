rule filter_pairs_hic:
	output:
		"resources/hic/files/pairs/{SRA_id}/{assembly_code}_{mapq_thr}.filtered.dedup.sorted.pairs"
	input:
		"resources/hic/files/pairs/{SRA_id}/{assembly_code}.dedup.sorted.pairs"
	log:
		"logs/filter_pairs_hic/{assembly_code}_{mapq_thr}_{SRA_id}.log"
	params:
		num_cores_hic=num_cores_hic
	conda:
		"../../envs/hic_processing_env.yml"
	shell:
		"""
		pairtools select "mapq1>={wildcards.mapq_thr} and mapq2>={wildcards.mapq_thr}" \
		--nproc-in {params.num_cores_hic} \
		--nproc-out {params.num_cores_hic} \
		{input} \
		-o {output} 2>{log}
		"""
