rule annotate_backbone:
	output:
		backbone_full="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}.filtered.sorted.annotated.bed",
		backbone_mask="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}.mask.filtered.sorted.annotated.bed"
	input:
		backbone_full="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.filtered.sorted.bed",
		backbone_mask="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.mask.filtered.sorted.bed",
		annotation="resources/reference_data/gencode/PC/v{gencode_version}_gene.sorted.bed"
	conda:
		"../../envs/bedtools.yml"
	log:
		"logs/annotate_backbone/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}.log"
	shell:
		"""
		bedtools intersect -a {input.backbone_full} -b {input.annotation} -wa -wb | \
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
		bedtools sort > {output.backbone_full} 2>> {log}

		bedtools intersect -a {input.backbone_mask} -b {input.annotation} -wa -wb | \
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
		bedtools sort > {output.backbone_mask} 2>> {log}
		"""
