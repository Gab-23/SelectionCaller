rule process_haploinsufficiency:
    output:
        "data/haploinsufficiency/{assembly_code}_{window_mode}_{window_size}_{step}_{haplo_thr}"
    input:
        backbone = "resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.filtered.sorted.annotated.bed",
        haplo = haploinsufficiency_location
    log:
        "logs/process_haploinsufficiency/{assembly_code}_{window_mode}_{window_size}_{step}_{haplo_thr}.log"
    params:
        C = num_cores_R
    conda:
        "../envs/r_spdep.yml"
    shell:
        """
        Rscript workflow/scripts/process_haploinsufficiency.R \
        --input_path_haplo = {input.path_haplo} \
        --input_path_backbone = {input.path_backbone} \
        --output_path = {output} \
        --haplo_thr = {wildcards.haplo_thr} \
        --C = {params.C}
        """