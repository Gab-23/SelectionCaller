checkpoint get_mapping_EID_tt:
    output:
        "resources/epigenomes/mapping_EID_tt.tsv"
    log:
        "logs/get_mapping_EID_tt/get_mapping_EID_tt.tsv"
    shell:
        # NOTE (fix): homologated with get_reference_genome's download-robustness fix - piping
        # `wget -qO- URL | ...` directly hid wget's own errors (-q) and, since this is plain
        # (non-gzipped) TSV, a truncated download wouldn't crash - cut/sed/awk would just
        # silently process a partial export, with no retry and no visible sign anything was
        # wrong. Now downloads to a real temp file first (with retry/timeout/continue), then
        # processes that file through the same cut/sed/awk chain.
        """
        TMP_TSV=$(mktemp)
        wget --tries=5 --timeout=60 --continue -O "$TMP_TSV" \
        "https://docs.google.com/spreadsheets/d/1yikGx4MsO9Ei36b64yOy9Vb6oPC5IBGlFbYEt-N6gOM/export?format=tsv&gid=15" \
        2>{log}
        cat "$TMP_TSV" | \
        cut -f 2,4,6,15,17,18 | \
        sed '/^[[:space:]]*$/d' | \
        awk -F"\t" 'BEGIN{{OFS="\t"}} NR==1{{print "epigenome_id", "group_description","epigenome_mnemonic","standard_name","tissue_type","sample_type"}} NR>1{{print $1,$2,$3,$4,$5,$6}}' \
        > {output} 2>>{log}
        rm -f "$TMP_TSV"
        """
