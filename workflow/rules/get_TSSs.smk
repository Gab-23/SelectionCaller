rule get_TSSs:
    output:
	    "resources/reference_data/gencode/{assembly_code}_v{gencode_version}_{tss_subset}.bed"
    input:
        "resources/reference_data/gencode/gencode.v{gencode_version}.annotation.gff3.gz"
    conda:
        "../envs/r_backbone.yml"
    log:
        "logs/get_TSSs/{assembly_code}_v{gencode_version}_{tss_subset}.log"
    shell:
        """
         Rscript workflow/scripts/get_TSSs.R \
        --input_path={input} \
        --output_path={output} 2> {log}
        """