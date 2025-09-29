rule get_epigenome_files:
    output:
        "resources/roadmap_epigenomes_download_log_{tissue_type}_{liftover_str}_{file_format}"
    input:
	    "resources/roadmap_epigenomes_index_{tissue_type}_{liftover_str}_{file_format}"  
    shell:
        r"""
        xargs -I{{}} sh -c '
        filename=$(echo "{{}}" | awk -v name="$(basename "{{}}")" '\''{{gsub(/\.gz$/,"", name); print name}}'\'');
        wget -qO- "{{}}" | gunzip > "data/$filename";
        echo "data/$filename" >> {output}
        ' < {input}
        """