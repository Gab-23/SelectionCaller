rule get_epigenome_files:
    output:
        "resources/epigenomes/checks/downloads/{tissue_type}_{sample_type}_{liftover_str}_{file_format}"
    input:
	    "resources/epigenomes/checks/paths/{tissue_type}_{sample_type}_{liftover_str}_{file_format}"
    log:
        "logs/get_epigenome_files/{tissue_type}_{sample_type}_{liftover_str}_{file_format}.log"
    shell:
        r"""
        xargs -I{{}} sh -c '
        filename=$(echo "{{}}" | awk -v name="$(basename "{{}}")" '\''{{gsub(/\.gz$/,"", name); print name}}'\'');
        wget -qO- "{{}}" | gunzip > "resources/epigenomes/roadmap/$filename";
        echo "resources/epigenomes/roadmap/$filename" >> {output}
        ' < {input} 2>{log}
        """