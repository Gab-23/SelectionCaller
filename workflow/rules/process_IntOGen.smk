rule process_IntOGen:
    output:
	    "resources/IntOGen/processed/{assembly_code}_v{gencode_version}_{tumor_type_IntOGen}.bed"
    input:
        path_IntOGen=config["IntOGen_locations"][0],
        path_protein_coding_transcript="resources/reference_data/gencode/PC/v{gencode_version}_transcript.bed",
    params:
        cohort=cohort_IntOGen
    conda:
        "../envs/r_roadmap.yml"
    log:
        "logs/process_IntOGen/{assembly_code}_v{gencode_version}_{tumor_type_IntOGen}.log"
    shell:
        """
        Rscript workflow/scripts/process_IntOGen.R \
        --input_path_IntOGen={input.path_IntOGen} \
        --input_path_protein_coding_transcript={input.path_protein_coding_transcript} \
        --tumor_type_IntOGen={wildcards.tumor_type_IntOGen} \
        --cohort={params.cohort} \
        --output_path={output} 2>{log}
        """