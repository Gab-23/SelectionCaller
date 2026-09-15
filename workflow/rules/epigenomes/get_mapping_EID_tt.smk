checkpoint get_mapping_EID_tt:
    output:
        "resources/epigenomes/mapping_EID_tt.tsv"
    input:
        # NOTE (fix): this used to download the sheet live via wget on every run, but the
        # cluster's outbound firewall blocks googleusercontent.com (the CDN host Google
        # Sheets' export redirects to) - confirmed directly with two independent wget attempts,
        # same host, same "Connection timed out" / "Network is unreachable" failure both
        # times, ruling out a one-off transient issue. No amount of retry/timeout fixes a
        # firewall block, so this is now a static, manually-downloaded resource instead (the
        # mapping table changes rarely) - see README for how to refresh it if it ever needs
        # updating.
        "resources/epigenomes/mapping_EID_tt_raw.tsv"
    log:
        "logs/get_mapping_EID_tt/get_mapping_EID_tt.tsv"
    shell:
        """
        cat {input} | \
        cut -f 2,4,6,15,17,18 | \
        sed '/^[[:space:]]*$/d' | \
        awk -F"\t" 'BEGIN{{OFS="\t"}} NR==1{{print "epigenome_id", "group_description","epigenome_mnemonic","standard_name","tissue_type","sample_type"}} NR>1{{print $1,$2,$3,$4,$5,$6}}' \
        > {output} 2>{log}
        """