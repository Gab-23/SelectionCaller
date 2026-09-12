rule get_epigenome_files:
    output:
        "resources/epigenomes/checks/downloads/{tissue_type_Roadmap}_{sample_type}_{liftover_str}_{file_format}"
    input:
        "resources/epigenomes/checks/paths/{tissue_type_Roadmap}_{sample_type}_{liftover_str}_{file_format}"
    log:
        "logs/get_epigenome_files/{tissue_type_Roadmap}_{sample_type}_{liftover_str}_{file_format}.log"
    shell:
        # NOTE (fix): homologated with get_reference_genome's download-robustness fix - each
        # iteration used to pipe `wget -qO- "{}" | gunzip` directly, hiding wget's own errors
        # (-q) and silently writing a truncated/corrupt file if any one of these (potentially
        # many) downloads was incomplete, with no retry. Now each URL is downloaded to a real
        # temp file first (with retry/timeout/continue), then decompressed as a separate step.
        r"""
        xargs -I{{}} sh -c '
        filename=$(echo "{{}}" | awk -v name="$(basename "{{}}")" '\''{{gsub(/\.gz$/,"", name); print name}}'\'');
        TMP_GZ=$(mktemp --suffix=.gz);
        wget --tries=5 --timeout=60 --continue -O "$TMP_GZ" "{{}}";
        gunzip -c "$TMP_GZ" > "resources/epigenomes/roadmap/$filename";
        rm -f "$TMP_GZ";
        echo "resources/epigenomes/roadmap/$filename" >> {output}
        ' < {input} 2>{log}
        """
