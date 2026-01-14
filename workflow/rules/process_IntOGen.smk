rule process_IntOGen:
    output:
	    "resources/IntOGen/processed/{assembly_code}_v{gencode_version}_{tss_subset}_{tumor_type_IntOGen}.bed"
    input:
        path_IntOGen=config["IntOGen_locations"][0],
        path_TSSs="resources/reference_data/gencode/{assembly_code}_v{gencode_version}_{tss_subset}.sorted.bed",
    params:
        cohort=config["cohorts"][0]
    conda:
        "../envs/r_roadmap.yml"
    log:
        "logs/process_IntOGen/{assembly_code}_v{gencode_version}_{tss_subset}_{tumor_type_IntOGen}.log"
    shell:
        """
        Rscript workflow/scripts/process_IntOGen.R \
        --input_path_IntOGen={input.path_IntOGen} \
        --input_path_TSSs={input.path_TSSs} \
        --tumor_type_IntOGen={wildcards.tumor_type_IntOGen} \
        --cohort={params.cohort} \
        --output_path={output} 
        """