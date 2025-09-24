rule get_epigenome_paths:
    output:
	    "resources/roadmap_epigenomes_index_{tumor_type}_{is_liftover}_{file_format}"  
    params:
        input_url = "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/imputed12marks/jointModel/final/",
        ids = lambda wildcards: "|".join(epigenome_ids),
        conditional_grep = lambda wildcards: 'grep -E hg38lift' if wildcards.is_liftover else ''
    shell:
        """
        wget -qO- {params.input_url} | \
        grep -oE '"E[0-9]+[^"]*"' | \
        awk -v prefix="{params.input_url}" '{{gsub(/"/, "", $1); print prefix $1}}' | \
        grep -E "{params.ids}" | \
        grep -E "{wildcards.file_format}.bed.gz" | \
        {params.conditional_grep} > {output}
        """