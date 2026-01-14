rule get_epigenome_paths:
    output:
	    "resources/epigenomes/checks/paths/{tissue_type_Roadmap}_{sample_type}_{liftover_str}_{file_format}"
    input:
        mapping_file=lambda wildcards: checkpoints.get_mapping_EID_tt.get().output[0]
    params:
        input_url = "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/imputed12marks/jointModel/final/",
        ids = lambda wildcards: "|".join(get_epigenome_ids(wildcards, tissue_type_Roadmap)),
        conditional_grep = lambda wildcards: '| grep -E hg38lift' if wildcards.liftover_str == 'hg38lift' else '| grep -vE hg38lift'
    log:
        "logs/get_epigenome_paths/{tissue_type_Roadmap}_{sample_type}_{liftover_str}_{file_format}.log"
    shell:
        """
        wget -qO- {params.input_url} | \
        grep -oE '"E[0-9]+[^"]*"' | \
        awk -v prefix="{params.input_url}" '{{gsub(/"/, "", $1); print prefix $1}}' | \
        grep -E "{params.ids}" | \
        grep -E "{wildcards.file_format}.bed.gz" \
        {params.conditional_grep} > {output} 2> {log}
        """