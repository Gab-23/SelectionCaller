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
        # NOTE (fix): homologated with get_reference_genome's download-robustness fix - piping
        # `wget -qO- URL | grep | ...` directly hid wget's own errors (-q) and, since this is
        # plain (non-gzipped) HTML, a truncated download wouldn't even crash downstream -
        # grep/awk would just silently process a partial page, with no retry and no visible
        # sign anything was wrong. Now downloads to a real temp file first (with retry/
        # timeout/continue), then processes that file through the same grep/awk chain.
        """
        TMP_HTML=$(mktemp)
        wget --tries=5 --timeout=60 --continue -O "$TMP_HTML" {params.input_url} 2> {log}
        cat "$TMP_HTML" | \
        grep -oE '"E[0-9]+[^"]*"' | \
        awk -v prefix="{params.input_url}" '{{gsub(/"/, "", $1); print prefix $1}}' | \
        grep -E "{params.ids}" | \
        grep -E "{wildcards.file_format}.bed.gz" \
        {params.conditional_grep} > {output} 2>> {log}
        rm -f "$TMP_HTML"
        """