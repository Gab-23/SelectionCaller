rule filter_manifest:
    output:
        "resources/GDC/CNAs/{tumor_type_GDC}_{gender_GDC}_manifest"
    input:
        path_manifest=config["GDC_manifest_locations"][0],
        path_files="resources/GDC/CNAs/{tumor_type_GDC}_{gender_GDC}_files"
    log:
        "logs/filter_manifest/{tumor_type_GDC}_{gender_GDC}.log"
    shell:
        """
        head -n 1 {input.path_manifest} > {output} 2>>{log}
        grep -Fw -f {input.path_files} {input.path_manifest} >> {output} 2>>{log}
        """