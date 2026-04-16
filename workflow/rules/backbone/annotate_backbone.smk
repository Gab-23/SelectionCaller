rule annotate_backbone:
	output:
		"resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.filtered.sorted.annotated.bed"
	input:
		backbone="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.filtered.sorted.bed",
		annotation="resources/reference_data/gencode/PC/v{gencode_version}_gene.sorted.bed"
	conda:
		"../envs/bedtools.yml"
	log:
		"logs/annotate_backbone/{assembly_code}_{window_mode}_{window_size}_{step}.log"
	shell:
		"""
		bedtools intersect -a {input.backbone} -b {input.annotation} -wa -wb | \
		cut -f 1,2,3,4,9 | \
		awk '
		{{
			key = $4
			if (genes[key] == "")
				genes[key] = $5
			else
				genes[key] = genes[key] "," $5
			start[key] = $2
			end[key] = $3
			chr[key] = $1
		}}
		END{{
			for (k in genes)
				print chr[k] "\t" start[k] "\t" end[k] "\t" k "\t" genes[k]
		}}
		' | \
		bedtools sort > {output} 2> {log}
		"""
