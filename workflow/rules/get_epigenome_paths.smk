rule get_epigenome_paths:
    input:
        mapping_file=lambda wildcards: checkpoints.get_mapping_eID_tt.get().output[0]
    output:
	    "resources/roadmap_epigenomes_index_{tissue_type}_{sample_type}_{liftover_str}_{file_format}"  
    params:
        input_url = "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/imputed12marks/jointModel/final/",
        ids = lambda wildcards: "|".join(get_epigenome_ids(wildcards)),
        conditional_grep = lambda wildcards: 'grep -E hg38lift' if wildcards.liftover_str else ''
    shell:
        """
        wget -qO- {params.input_url} | \
        grep -oE '"E[0-9]+[^"]*"' | \
        awk -v prefix="{params.input_url}" '{{gsub(/"/, "", $1); print prefix $1}}' | \
        grep -E "{params.ids}" | \
        grep -E "{wildcards.file_format}.bed.gz" | \
        {params.conditional_grep} > {output}
        """