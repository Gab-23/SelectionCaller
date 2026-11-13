rule compute_cds_overlap:
    # bp di sequenza codificante per bin: CDS dei geni protein_coding di GENCODE, fuse tra
    # trascritti (ogni bp contato una volta). Serve come esposizione per le mutazioni
    # codificanti (filtro dei bin e densita' per kb), non come predittore.
    output:
        "data/num_overlaps/CDS/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}"
    input:
        path_gtf="resources/reference_data/gencode/gencode.v{gencode_version}.annotation.gtf.gz",
        path_backbone="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.mask.filtered.sorted.bed"
    conda:
        "../../envs/r_roadmap.yml"
    log:
        "logs/compute_cds_overlap/{assembly_code}_{window_mode}_{window_size}_{step}_v{gencode_version}.log"
    shell:
        """
        Rscript workflow/scripts/compute_cds_overlap.R \
        --input_path_gtf={input.path_gtf} \
        --input_path_backbone={input.path_backbone} \
        --output_path={output} 2>{log}
        """